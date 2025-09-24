rm(list = ls())

library("parallel")
library("data.table")
library("dplyr")

#chromosomes = c(seq(1,22))
chromosomes = c(11)


options(width = 180)

#Define the file names.
storage_directory = "/home/sharedFolder/projects/geneticConfounding/conference/data/SickleCellDosage/"
#samples_to_include = paste0(storage_directory, "T2DCohort.tsv")
#ukbb_genotype_source_dir = "/home/sharedFolder/referenceData/ukb/genotypes/"
    

extractCohort = function(chromosome = NULL, storage_directory, samples_to_include, ukbb_genotype_source_dir) {
    print(paste("Creating extracted CHR file for chromosome", chromosome))
    ukbb_genotype_chr_file_prefix = paste(ukbb_genotype_source_dir,
                                   'ukb22418_c',
                                   chromosome,
                                   '_b0_v2',
                                    sep = "")
    
    ukbb_genotype_extracted_file_path = paste(storage_directory,
                                             "extractedChr",
                                             chromosome,
                                             sep = "")
    
    ukbb_genotype_ld_prune_variants = paste(storage_directory,
                                             "extractedChr",
                                             chromosome,
                                             ".LDSNPs",
                                             sep = "")
    
    ukbb_genotype_ld_pruned_file_path = paste(storage_directory,
                                             "extractedChr",
                                             chromosome,
                                             ".LDPruned",
                                             sep = "")
    
    #Extract files by chrs.
    command = paste('plink',
        '--bfile',
        ukbb_genotype_chr_file_prefix,        
        '--keep-allele-order',
        '--allow-no-sex',
        '--keep-fam',
        samples_to_include,     
        '--recode12',
        'transpose',
        '--out', ukbb_genotype_extracted_file_path)
    print(command)
    system(command)
    
    #Get LD SNPs.
    command = paste('plink',
        '--tfile',
        ukbb_genotype_extracted_file_path,        
        '--keep-allele-order',
        '--allow-no-sex',
        '--indep-pairwise 50 10 0.1',
        '--out', ukbb_genotype_ld_prune_variants)
    print(command)
    system(command)
    
    #Get LD SNPs.
    command = paste('plink',
        '--tfile',
        ukbb_genotype_extracted_file_path,        
        '--keep-allele-order',
        '--allow-no-sex',
        '--extract',
        paste0(ukbb_genotype_ld_prune_variants, ".prune.in"),
        '--recode12',
        'transpose',
        '--out', ukbb_genotype_ld_pruned_file_path)
    print(command)
    system(command)
    
    
}

#Don't run this more than once.
#mclapply(chromosomes, function(i) extractCohort(i, storage_directory, samples_to_include, ukbb_genotype_source_dir), mc.cores = 22)

print("All cohort genotype information is extracted.")
print("Creating dosage files now...")
         
###Get the dosage file now.
createDosageFile = function(chromosome = NULL, storage_directory) {
    #tped_file = paste0(storage_directory, "extractedChr", chromosome, ".LDPruned.tped")
    #tfam_file = paste0(storage_directory, "extractedChr", chromosome, ".LDPruned.tfam")
    
    #dosage_file = paste0(storage_directory, "extractedChr", chromosome, ".LDPruned.dosage.tsv")
    
    tped_file = paste0(storage_directory, "extractedBlackWhiteCohortForSickleCell", ".tped")
    tfam_file = paste0(storage_directory, "extractedBlackWhiteCohortForSickleCell", ".tfam")
    
    dosage_file = paste0(storage_directory, "extractedBlackWhiteCohortForSickleCell", ".dosage.tsv")
    
    #The following AWK command will sum every two haplotypes of a sample to get the dosage
    calculate_dosage_command = paste(
        "cat", 
        tped_file,
        "|",
        'awk -F " "',
        '\'{ printf($1"\t"$2"\t"$3"\t"$4"\t");',
        'for (i=5;i<=NF;i+=2) { printf($i+$(i+1)"\t") } print "" }\' ',
        '>',
        dosage_file)
    
    print(calculate_dosage_command)
    system(calculate_dosage_command)
}
         
mclapply(chromosomes, function(i) createDosageFile(i, storage_directory), mc.cores = 22)