rm(list=ls())
set.seed(13)

library("parallel")
library("reshape2")
library("data.table")
library("ggplot2")
library('dplyr')
library('mediation')
library('tidyverse')

options(scipen=999)

setwd("/home/sharedFolder/projects/geneticConfounding/scripts/")

chromosomes = c(4, 17, 3, 13, 20, 18, 1, 22, 11, 6, 2, 14, 16, 9, 10)

performMediationAnalysis = function(chromosome = NULL, merged_dataset, storage_dir) {
    print(paste("Working on chr: ",
               chromosome))
    dosage_file_path = paste0(storage_dir,
                            "extractedChr",
                            chromosome,
                            ".LDPruned.dosage.fixedHeader.tsv")
    dosage_data = as.data.frame(fread(dosage_file_path))
    
    #Put dbSNP IDs as indices.
    row.names(dosage_data) = dosage_data$dbSNP
    
    #Bifurcate the matrix.
    #Create the output dosage matrix with mediation effects set to 0 for now.
    dosage_snp_data = dosage_data[, 1:4]
    dosage_snp_data$AverageMediation = 0
    dosage_snp_data$AverageMediation.P = 0
    dosage_snp_data$AverageMediation.LCI = 0
    dosage_snp_data$AverageMediation.UCI = 0
    
    #Take the dosage values out for modeling.
    dosage_values = dosage_data[, 5:length(colnames(dosage_data))]
    
    #Replace all dosages with NA if they are 0. They are 0 only when they had mising genotype.
    dosage_values[dosage_values == 0] <- NA
    dosage_values = dosage_values - 2
    
    #Transpose the dosage matrix.
    dosage_values_transposed = as.data.frame(t(dosage_values))
    dosage_values_transposed$SampleIDs = row.names(dosage_values_transposed)
    
    merged_dataset = merge(merged_dataset, dosage_values_transposed, by.x = "eid", by.y = "SampleIDs")
    
    #Modify data variables.
    merged_dataset$CaseStatus = as.factor(merged_dataset$CaseStatus)
    
    #Remove NA entries for SED.
    merged_dataset = merged_dataset[!is.na(merged_dataset$SED), ]
    
    #Standardize and normalize the SED and PCs.
    merged_dataset$StandardizedSED = (merged_dataset$SED - mean(merged_dataset$SED)) / sd(merged_dataset$SED)
    
    merged_dataset$EthnicGroup = factor(merged_dataset$EthnicGroup, levels = c("White", "Black"))
    
    ############################################################################################################################
    
    mediation.totaleffect = glm(CaseStatus ~ SED + Age + Sex, 
                                family=binomial(link='logit'), data=merged_dataset)
    summary(mediation.totaleffect)
    
    print(paste("Running mediation on chr", chromosome))
    pb = txtProgressBar(min = 8, max = length(colnames(dosage_data)), initial = 0) 
    
    number_of_snps = length(colnames(merged_dataset)) - 1
    #print(number_of_snps)
    ###Mediation Analysis
    for (i in c(8:number_of_snps)){
        this_dataset = merged_dataset[!is.na(merged_dataset[[i]]),]
        
        #The mediator model
        #print(colnames(this_dataset)[i])
        mediation.mediator.Dosage.formula = paste(colnames(this_dataset)[i],
                                                  "~",
                                                  "SED + Age + Sex")
        
        mediation.mediator.Dosage = lm(mediation.mediator.Dosage.formula, data=this_dataset)
        
        #Adjusted model.
        mediation.dv.formula = paste("CaseStatus",
                                     "~",
                                     "SED + ",
                                     colnames(this_dataset)[i],
                                     " + Age + Sex")
        mediation.dv = glm(mediation.dv.formula, 
                           family=binomial(link='logit'), data=this_dataset)
        #print(summary(mediation.dv))
        
        #Check for linear dependence of predictior variables in GLM; if found, don't run the mediation model.
        #Check this out for more: https://stats.stackexchange.com/questions/25839/logistic-regression-in-r-returning-na-values
        if (anyNA(mediation.dv$coefficients)) {
            dosage_snp_data[colnames(this_dataset)[i],][5:8] = c(-9,-9,-9,-9)
            next
        }
        
        #Mediation
        results = mediate(mediation.mediator.Dosage, mediation.dv, 
                          treat='SED', 
                          mediator=paste(colnames(this_dataset)[i]),
                          covariates = c("Age", "Sex"),
                          sims=10)
        
        dosage_snp_data[colnames(this_dataset)[i],][5:8] = c(results$n.avg, 
                                                             results$n.avg.p,
                                                             results$n.avg.ci[[1]],
                                                             results$n.avg.ci[[2]])
        
        #print(summary(results))
        #Set up progress bar report.
        setTxtProgressBar(pb,i)

    }
    #Close progress bar.
    close(pb)
    
    print(paste("Writing mediation results on chr", chromosome))
    write.table(dosage_snp_data, 
                file = paste0(storage_dir,
                              "extractedChr",
                              chromosome,
                              ".mediation.tsv"),
                sep = "\t",
                row.names = F,
                quote = F)
}

print("Reading demographic and SED files now...")
storage_dir = "../data/UKBBSelectCohort/T2D/"

#Read the demographic data.
demographics = as.data.frame(fread("../data/ukb_demographics.tsv"))

sed = as.data.frame(fread("../data/UKBBSed.tsv"))
case_control_cohort = as.data.frame(fread("../data/PheCodeFiles/250.2_Type_2_diabetes.txt", col.names = c("EID", "CaseStatus")))


merged_dataset = merge(demographics, sed, by.x = "eid", by.y = "Individual")
merged_dataset = merge(merged_dataset, case_control_cohort, by.x = "eid", by.y = "EID")

mclapply(chromosomes, function(i) performMediationAnalysis(i, 
                                                merged_dataset, storage_dir), mc.cores = 22)
