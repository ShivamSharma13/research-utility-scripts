rm(list=ls())
set.seed(13)

library("data.table")
library("dplyr")
library("parallel")
library("ggplot2")
library("tidyr")

options(repr.plot.width=24, repr.plot.height=20)

#Read the phecode map file which we used to create the phecode cohorts.
phecode_icd10_map = as.data.frame(fread("Data/PhecodeFiles/ReferenceFiles/Phecode_map_v1_2_icd10cm_beta.csv"))

#Get a phecodes we have the data for. They are 1,755.
phecode_icd10_map = phecode_icd10_map %>% select(phecode) %>%distinct(phecode)

head(phecode_icd10_map)
nrow(phecode_icd10_map)

###Read the phecodes definition file.
phecode_definitions = as.data.frame(fread("Data/PhecodeFiles/ReferenceFiles/phecode_definitions1.2.csv"))

head(phecode_definitions)
nrow(phecode_definitions)

#Merge the definitions witht the ICD10CM map.
phecode_information = merge(phecode_icd10_map, phecode_definitions, by.x = "phecode", by.y = "phecode")

#Get only the relevant columns from Phecodes.
phecode_information = phecode_information %>% select(phecode, phenotype, category, sex)

head(phecode_information)
nrow(phecode_information)

#Read the Rye estimates file.
rye_estimates = as.data.frame(fread("Data/Rye/extractedChrAllPruned.ContinentalEstimates.6.Q"))
rye_estimates = rye_estimates %>% select(person_id, SelfReportedRaceEthnicity, European, African, NativeAmerican)

#Also bring the demographics and merge with ancestry df.
demographics = as.data.frame(fread("Data/Metadata//DemographicData.tsv"))

demographics = demographics %>% select(person_id, sex_at_birth, date_of_birth)

#Calculate age of each samples.
demographics$Age = as.numeric(difftime(Sys.Date(),demographics$date_of_birth, units = "weeks"))/52.25

#Select only Males & Females for this analysis.
demographics = demographics %>% filter(sex_at_birth %in% c("Male", "Female")) %>% select(person_id, sex_at_birth, Age)
names(demographics) = c("person_id", "Sex", "Age")

nrow(rye_estimates)
nrow(demographics)

#Merge demographics and ancestry estimates.
metadata = merge(rye_estimates, demographics, by = "person_id")

head(metadata)
nrow(metadata)

#Get the number of samples with Genotype & EHR information.
sample_phecode_file = as.data.frame(fread("Data/PhecodeFiles/NewCohorts/cohort.585.2.tsv"))
head(sample_phecode_file)

#Merge them.
ancestry_phecode_merged_df = merge(metadata, sample_phecode_file, by.x = "person_id", by.y = "PersonID")
nrow(ancestry_phecode_merged_df)

phecode_cohort_directory = "/home/jupyter/workspaces/geneticancestry/Data/PhecodeFiles/NewCohorts/"


getNumberOfCasesForPhecodes = function (phecode){
    cohort_file_name = paste("cohort",phecode,"tsv",sep = ".")
    cohort_file_path = paste(phecode_cohort_directory, cohort_file_name, sep="/")
    
    #Read the phecode cohort file.
    phecode_cohort = as.data.frame(fread(cohort_file_path))
    
    #Merge with Metadata and ancestry data.
    complete_cohort = merge(metadata, phecode_cohort, by.x = "person_id", by.y = "PersonID")
    return(paste(nrow(complete_cohort[complete_cohort$Status == "Case",]),
                nrow(complete_cohort[complete_cohort$Status == "Control",]),
                nrow(complete_cohort[complete_cohort$Status == "Excluded",]),
                nrow(complete_cohort), sep = "-"))
}

#Get cohort counts for each phecode.
phecode_information = phecode_information %>% 
    rowwise() %>% 
    mutate(CohortNumbers = getNumberOfCasesForPhecodes(phecode))

head(phecode_information)

#Extract the cohort numbers for each Phecode.
phecode_information = phecode_information %>%
  separate(CohortNumbers, c("Case", "Control", "Excluded", "Total"), "-")

phecode_information$Prevalence = as.numeric(phecode_information$Case)/ as.numeric(phecode_information$Total) * 100
head(phecode_information)

#Write the Phecode information file.
write.table(file = "Data/PhecodeFiles/PhecodeSummarizedInfo.tsv", x = phecode_information, row.names = F, quote = F, sep = "\t")

phecode_information = as.data.frame(fread("Data/PhecodeFiles/PhecodeSummarizedInfo.tsv"))

ggplot(phecode_information, aes(x=Prevalence)) + 
    geom_histogram(color="black", fill="gray") +
    xlab("Prevalence Percentage") +
    xlim(0,100) + theme_classic(base_size = 32)

ggplot(phecode_information, aes(x=Prevalence)) + 
    geom_histogram(color="black", fill="gray") +
    xlab("Prevalence Percentage") +
    xlim(0,10) + theme_classic(base_size = 32)

ggplot(phecode_information, aes(x=Prevalence)) + 
    geom_histogram(color="black", fill="gray") +
    xlab("Prevalence Percentage") +
    xlim(0,1) + theme_classic(base_size = 32)

###Break the prevalence into buckets for further analysis.
phecode_information$PrevalenceLimits = cut(phecode_information$Prevalence, breaks=c(-0.001, 0, 0.01, 0.02, 0.05, 
                                                                                    0.1, 0.5, 1, 
                                                                                    5, 10, 15, 20, 25, 
                                                                                    30, 35), 
                                           labels = c("0", "0-0.01", "0.01-0.02", "0.02-0.05",
                                                      "0.05-0.1", "0.1-0.5", "0.5-1",
                                                      "1-5", "5-10", "10-15", "15-20", "20-25",
                                                      "25-30", "30-35"),
                                            include.lowest = TRUE)

prevalence_limits_phcodes = count(phecode_information, PrevalenceLimits)

#Plot the prevalence bounds.
ggplot(data=prevalence_limits_phcodes, aes(x=PrevalenceLimits, y=n)) +
    geom_bar(stat="identity", fill="violet", color = "black") +
    geom_text(aes(label=n), vjust=-1.6, 
            color="black", size=8)+
    xlab("Prevalence percentage ranges") + ylab("Number of phecodes in the range") +
    theme_classic(base_size = 32) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(phecode_information[as.numeric(phecode_information$Case) < 100, ], aes(x=Case)) + 
    geom_histogram(color="black", fill="violet", bins = 10) +
   theme_classic(base_size = 32)


##We will ignore phecodes which have a prevelance of less than 0.01%
print(paste("Before filtering using prevelence", nrow(phecode_information)))

#Filter based on prevalence.
phecode_information = phecode_information[phecode_information$Case > 100,]
print(paste("After filtering using prevelence", nrow(phecode_information)))

summary_logits_columns = c("Phecode","AncestryEstimate","AncestryStdErr","AncestryZvalue","AncestryPvalue",
            "AgeEstimate","AgeStdErr","AgeZvalue","AgePvalue",
            "SexMaleEstimate","SexMaleStdErr","SexMaleZvalue","SexMalePvalue")

all_logits_africans = data.frame(matrix(nrow=1, ncol = length(summary_logits_columns)))
all_logits_europeans = data.frame(matrix(nrow=1, ncol = length(summary_logits_columns)))
all_logits_nativeamericans = data.frame(matrix(nrow=1, ncol = length(summary_logits_columns)))


colnames(all_logits_africans) = summary_logits_columns
colnames(all_logits_europeans) = summary_logits_columns
colnames(all_logits_nativeamericans) = summary_logits_columns

phecode_cohort_directory = "/home/jupyter/workspaces/geneticancestry/Data/PhecodeFiles/NewCohorts/"

regressDiseaseByAncestry = function(phecode, phecode_sex){
    print("-----------------------------------------------------------------------------------------")
    
    sex_specific_flag = FALSE
    cohort_file_name = paste("cohort", phecode, "tsv" ,sep = ".")
    cohort_file_path = paste(phecode_cohort_directory, cohort_file_name, sep="/")
    
    #Read the phecode cohort file.
    phecode_cohort = as.data.frame(fread(cohort_file_path))
    
    #Merge with metadata. This will leave us with samples who have ancestry data.
    phecode_cohort = merge(phecode_cohort, metadata, by.x = "PersonID", by.y = "person_id")
    
    #Get only the cases and controls.
    phecode_cohort = phecode_cohort[phecode_cohort$Status != "Excluded",]
    phecode_cohort$DiseaseStatus = ifelse(phecode_cohort$Status == "Case", 1, 0)
    
    print(paste0("Starting work on phecode: ", phecode))
    print(paste0("Phecode sex requirement is: ", phecode_sex))
    print(nrow(phecode_cohort))    
    #print(head(phecode_cohort))
    #print(phecode_cohort %>% count(Sex))
    
    #Filter if sex requirement is different.
    if (phecode_sex != "Both" & phecode_sex != "") {
        #Phecode is related to only Male or Female.
        print(paste0("Sex specific Phecode detected: ", phecode))
        print(paste0("Original cohort size: ", nrow(phecode_cohort)))
        
        #Get only the sex where this condition exists.
        phecode_cohort = phecode_cohort[phecode_cohort$Sex == phecode_sex,]
        print(paste0("New cohort size: ", nrow(phecode_cohort)))
        
        sex_specific_flag = TRUE
    }
        
    ###############################################################################################
    ###############################################################################################
    ###Modeling for Africans.
    if (sex_specific_flag == FALSE) {
        print("Running regular logit...")
        this_logit_african <- glm(DiseaseStatus ~ as.numeric(African/100) + Age + Sex, 
                             data = phecode_cohort, 
                             family = "binomial")
        #print(this_logit_african)
        #print(summary(this_logit_african))
        this_logit_summary_african_object = summary(this_logit_african)
        
        this_logit_summary_african = c(as.character(phecode),
                coef(this_logit_summary_african_object)["as.numeric(African/100)",1],
                    coef(this_logit_summary_african_object)["as.numeric(African/100)",2],
                    coef(this_logit_summary_african_object)["as.numeric(African/100)",3],
                    coef(this_logit_summary_african_object)["as.numeric(African/100)",4],
                coef(this_logit_summary_african_object)["Age",1],
                    coef(this_logit_summary_african_object)["Age",2],
                    coef(this_logit_summary_african_object)["Age",3],
                    coef(this_logit_summary_african_object)["Age",4],
                coef(this_logit_summary_african_object)["SexMale",1],
                    coef(this_logit_summary_african_object)["SexMale",2],
                    coef(this_logit_summary_african_object)["SexMale",3],
                    coef(this_logit_summary_african_object)["SexMale",4])
    }
    
    else {
        print("Running sex specific logit...")
        this_logit_african <- glm(DiseaseStatus ~ as.numeric(African/100) + Age, 
                             data = phecode_cohort, 
                             family = "binomial")
        #print(this_logit_african)
        #print(summary(this_logit_african))
        this_logit_summary_african_object = summary(this_logit_african)

        this_logit_summary_african = c(as.character(phecode),
                coef(this_logit_summary_african_object)["as.numeric(African/100)",1],
                    coef(this_logit_summary_african_object)["as.numeric(African/100)",2],
                    coef(this_logit_summary_african_object)["as.numeric(African/100)",3],
                    coef(this_logit_summary_african_object)["as.numeric(African/100)",4],
                coef(this_logit_summary_african_object)["Age",1],
                    coef(this_logit_summary_african_object)["Age",2],
                    coef(this_logit_summary_african_object)["Age",3],
                    coef(this_logit_summary_african_object)["Age",4],
                NaN, NaN, NaN, NaN)
    }
    #print(this_logit_summary_african)

    #Combine the data with the overall results df.
    all_logits_africans = rbind(all_logits_africans, this_logit_summary_african)
    assign('all_logits_africans',all_logits_africans,envir=.GlobalEnv)
    
    ###############################################################################################
    ###############################################################################################
    ###Modeling for Europeans.
    
    if (sex_specific_flag == FALSE) {
        this_logit_european <- glm(DiseaseStatus ~ as.numeric(European/100) + Age + Sex, 
                             data = phecode_cohort, 
                             family = "binomial")
        #print(this_logit_european)
        #print(summary(this_logit_european))
        this_logit_summary_european_object = summary(this_logit_european)
        
        this_logit_summary_european = c(as.character(phecode),
                coef(this_logit_summary_european_object)["as.numeric(European/100)",1],
                    coef(this_logit_summary_european_object)["as.numeric(European/100)",2],
                    coef(this_logit_summary_european_object)["as.numeric(European/100)",3],
                    coef(this_logit_summary_european_object)["as.numeric(European/100)",4],
                coef(this_logit_summary_european_object)["Age",1],
                    coef(this_logit_summary_european_object)["Age",2],
                    coef(this_logit_summary_european_object)["Age",3],
                    coef(this_logit_summary_european_object)["Age",4],
                coef(this_logit_summary_european_object)["SexMale",1],
                    coef(this_logit_summary_european_object)["SexMale",2],
                    coef(this_logit_summary_european_object)["SexMale",3],
                    coef(this_logit_summary_european_object)["SexMale",4])
    }
    
    else {
        this_logit_european <- glm(DiseaseStatus ~ as.numeric(European/100) + Age, 
                             data = phecode_cohort, 
                             family = "binomial")
        #print(this_logit_european)
        #print(summary(this_logit_european))
        this_logit_summary_european_object = summary(this_logit_european)
        
        this_logit_summary_european = c(as.character(phecode),
                coef(this_logit_summary_european_object)["as.numeric(European/100)",1],
                    coef(this_logit_summary_european_object)["as.numeric(European/100)",2],
                    coef(this_logit_summary_european_object)["as.numeric(European/100)",3],
                    coef(this_logit_summary_european_object)["as.numeric(European/100)",4],
                coef(this_logit_summary_european_object)["Age",1],
                    coef(this_logit_summary_european_object)["Age",2],
                    coef(this_logit_summary_european_object)["Age",3],
                    coef(this_logit_summary_european_object)["Age",4],
                NaN, NaN, NaN, NaN)
        
    }
    
    #print(this_logit_summary_european)
    
    #Combine the data with the overall results df.
    all_logits_europeans = rbind(all_logits_europeans, this_logit_summary_european)
    assign('all_logits_europeans',all_logits_europeans,envir=.GlobalEnv)
    
    ###############################################################################################
    ###############################################################################################
    ###Modeling for Native americans
    if (sex_specific_flag == FALSE) {
        this_logit_nativeamerican <- glm(DiseaseStatus ~ as.numeric(NativeAmerican/100) + Age + Sex, 
                             data = phecode_cohort, 
                             family = "binomial")
        #print(this_logit_nativeamerican)
        #print(summary(this_logit_nativeamerican))
        this_logit_summary_nativeamerican_object = summary(this_logit_nativeamerican)
        
        this_logit_summary_nativeamerican = c(as.character(phecode),
                coef(this_logit_summary_nativeamerican_object)["as.numeric(NativeAmerican/100)",1],
                    coef(this_logit_summary_nativeamerican_object)["as.numeric(NativeAmerican/100)",2],
                    coef(this_logit_summary_nativeamerican_object)["as.numeric(NativeAmerican/100)",3],
                    coef(this_logit_summary_nativeamerican_object)["as.numeric(NativeAmerican/100)",4],
                coef(this_logit_summary_nativeamerican_object)["Age",1],
                    coef(this_logit_summary_nativeamerican_object)["Age",2],
                    coef(this_logit_summary_nativeamerican_object)["Age",3],
                    coef(this_logit_summary_nativeamerican_object)["Age",4],
                coef(this_logit_summary_nativeamerican_object)["SexMale",1],
                    coef(this_logit_summary_nativeamerican_object)["SexMale",2],
                    coef(this_logit_summary_nativeamerican_object)["SexMale",3],
                    coef(this_logit_summary_nativeamerican_object)["SexMale",4])
    }
    
    else {
        this_logit_nativeamerican <- glm(DiseaseStatus ~ as.numeric(NativeAmerican/100) + Age, 
                             data = phecode_cohort, 
                             family = "binomial")
        #print(this_logit_nativeamerican)
        #print(summary(this_logit_nativeamerican))
        this_logit_summary_nativeamerican_object = summary(this_logit_nativeamerican)
        
        this_logit_summary_nativeamerican = c(as.character(phecode),
                coef(this_logit_summary_nativeamerican_object)["as.numeric(NativeAmerican/100)",1],
                    coef(this_logit_summary_nativeamerican_object)["as.numeric(NativeAmerican/100)",2],
                    coef(this_logit_summary_nativeamerican_object)["as.numeric(NativeAmerican/100)",3],
                    coef(this_logit_summary_nativeamerican_object)["as.numeric(NativeAmerican/100)",4],
                coef(this_logit_summary_nativeamerican_object)["Age",1],
                    coef(this_logit_summary_nativeamerican_object)["Age",2],
                    coef(this_logit_summary_nativeamerican_object)["Age",3],
                    coef(this_logit_summary_nativeamerican_object)["Age",4],
                NaN, NaN, NaN, NaN)   
    }
    
    #print(this_logit_summary_nativeamerican)
    
    #Combine the data with the overall results df.
    all_logits_nativeamericans = rbind(all_logits_nativeamericans, this_logit_summary_nativeamerican)
    assign('all_logits_nativeamericans',all_logits_nativeamericans,envir=.GlobalEnv)
    
    ###############################################################################################
    
    
    return(0)
}

#Run the logistic model function for each row (each row has info about one Phecode)
phecode_information %>% 
    rowwise() %>% 
    mutate(Useless = regressDiseaseByAncestry(phecode, sex))

#Remove the dummy rows from the dfs.
head(all_logits_africans, 20)
head(all_logits_europeans, 20)
head(all_logits_nativeamericans, 20)

#Remove the dummy rows from the dfs.
nrow(all_logits_africans)
nrow(all_logits_europeans)
nrow(all_logits_nativeamericans)

#Write the dfs to a file.
write.table(file = "Data/PhecodeFiles/LogisticModels/NewCohort/African.tsv", x = all_logits_africans, 
            row.names = F, quote = F, sep = "\t")

write.table(file = "Data/PhecodeFiles/LogisticModels/NewCohort/European.tsv", x = all_logits_europeans, 
            row.names = F, quote = F, sep = "\t")

write.table(file = "Data/PhecodeFiles/LogisticModels/NewCohort/NativeAmerican.tsv", x = all_logits_nativeamericans, 
            row.names = F, quote = F, sep = "\t")

head(summary_logits, 20)

#Make a logit model for just T2D.
#Read the phecode cohort file.
phecode_cohort = as.data.frame(fread("/home/jupyter/workspaces/geneticancestry/Data/PhecodeFiles/Cohorts/cohort.250.2.tsv"))

#Merge with metadata.
phecode_cohort = merge(phecode_cohort, metadata, by.x = "PersonID", by.y = "person_id")

#Get only the cases and controls.
phecode_cohort = phecode_cohort[phecode_cohort$Status != "Excluded",]
phecode_cohort$DiseaseStatus = ifelse(phecode_cohort$Status == "Case", 1, 0)

#Run the logistic model
logit <- glm(DiseaseStatus ~ as.numeric(African) + Age + Sex, data = phecode_cohort, family = "binomial")

summary(logit)

#Get odds ratio.
oddsratio = exp(cbind(coef(logit), confint(logit)))  

print(oddsratio)

oddsratio["African",3]


































