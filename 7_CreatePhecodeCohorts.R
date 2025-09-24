rm(list=ls())
set.seed(13)

library("data.table")
library("dplyr")
library("parallel")
library('stringr')

options(repr.matrix.max.cols=50, repr.matrix.max.rows=100)

#Read the cohort file.
disease_cohort = as.data.frame(fread("Data/Metadata/AllDiseases.tsv"))
head(disease_cohort)

#Number of total diagnosis.
nrow(disease_cohort)

#Number of unique participants.
length(unique(disease_cohort$person_id))

#Diagnosis with only ICD10cm or ICD9CM info.
disease_cohort_with_icds = disease_cohort[(disease_cohort$source_vocabulary == "ICD10CM") | 
                                     (disease_cohort$source_vocabulary == "ICD9CM"),]

#Look at the number of samples for which we can get a Phecode information.
length(unique(disease_cohort_with_icds$person_id))
head(disease_cohort_with_icds)

#Get the date of assesment to implement the rule of two.
disease_cohort_with_icds$DateOfAssesment = format(disease_cohort_with_icds$condition_start_datetime, format = "%Y~%m~%d")

#Look at the number of samples for which we can get a Phecode information.
length(unique(disease_cohort_with_icds$person_id))
head(disease_cohort_with_icds)

#Let's look at number of distinct samples in each of source vocabulary.
disease_cohort %>% group_by(source_vocabulary) %>% summarize(n_distinct(person_id))

icd10cm = disease_cohort %>% filter(source_vocabulary == "ICD10CM") %>% select(person_id) %>% distinct(person_id)
icd9cm = disease_cohort %>% filter(source_vocabulary == "ICD9CM") %>% select(person_id) %>% distinct(person_id)
snomed = disease_cohort %>% filter(source_vocabulary == "SNOMED") %>% select(person_id) %>% distinct(person_id)

print(paste("Samples in icd10cm and not in icd9cm", nrow(setdiff(icd10cm, icd9cm))))
print(paste("Samples in icd10cm and not in snomed", nrow(setdiff(icd10cm, snomed))))

print(paste("Samples in icd9cm and not in icd10cm", nrow(setdiff(icd9cm, icd10cm))))
print(paste("Samples in snomed and not in icd10cm", nrow(setdiff(snomed, icd10cm))))

print("Sampling 1 million random code diagnosis.")
disease_cohort %>% sample_n(1000000) %>% group_by(source_vocabulary) %>% summarize(n_distinct(person_id))

print("Avergae date for all diagnosis by Source Codes")
disease_cohort %>% sample_n(1000000) %>% group_by(source_vocabulary) %>% summarize(MeanDate = mean(condition_start_datetime))

#Read the phecode map file.
phecode_icd10_map = as.data.frame(fread("Data/PhecodeFiles//ReferenceFiles/Phecode_map_v1_2_icd10cm_beta.csv"))
head(phecode_icd10_map)

#Read the phecode map file.
all_icd_code_to_phecode_map = as.data.frame(fread("Data/PhecodeFiles//ReferenceFiles/ICD-CMToPheCodeUnrolled.txt"))

#Get the phecode ICD9 map.
phecode_icd9_map = all_icd_code_to_phecode_map[all_icd_code_to_phecode_map$flag == "9",]
head(phecode_icd9_map)

#Get a list of ICD10 codes needed to be included for each PheCode.
phecode_map_icd_10_aggregated = phecode_icd10_map %>% 
                            group_by(phecode) %>% 
                            summarise(IncludeICD10s = list(icd10cm), ExcludePheCodeRange=unique(exclude_range))
head(phecode_map_icd_10_aggregated)

#Get the include ICD9 codes for each PheCodes.
phecode_icd9_map = phecode_icd9_map %>% 
                    dplyr::filter(flag == "9") %>%
                    group_by(phecode) %>% 
                    summarise(IncludeICD9s = list(ICD))

#Select only the Phecode and ICD 9 codes to include
phecode_map_icd_9_aggregated = phecode_icd9_map %>% 
                    select(phecode, IncludeICD9s)

head(phecode_map_icd_9_aggregated)

#Merge the two ICD maps and keep all rows from ICD10 maps.
phecode_map_aggregated = merge(phecode_map_icd_10_aggregated, phecode_map_icd_9_aggregated, 
                           by.x = "phecode", by.y = "phecode",
                           all.x = TRUE)

head(phecode_map_aggregated,4)

#Look at the number of rows in both merged & ICD10 maps.
print(nrow(phecode_map_aggregated))
print(nrow(phecode_map_icd_10_aggregated))

###Let's create the phecode files here.
getExclusionICD10s = function(phecode, include_icds, phecodes_exclude_ranges, include_icd_column_name){
    #Since phecode exclusion criteria can have multiple ranges, we assumer there are always more than one.
    phecodes_exclude_ranges = as.list(strsplit(phecodes_exclude_ranges, ",")[[1]])
    
    #Global list of exclusion ICD10 or ICD9s codes to be returned.
    exclude_icds = c()
    
    #Get ranges of Phecodes that need to be excluded.
    for (this_excluded_range in phecodes_exclude_ranges){
        #We extract the upper and lower bounds for this exclusion range.
        this_excluded_range_bounds = as.list(strsplit(this_excluded_range, "-")[[1]])
        
        #Find all the Phecodes that lie between this exclusion range.
        excluded_phecodes = phecode_map_aggregated %>% 
            filter((as.numeric(phecode) >= as.numeric(this_excluded_range_bounds[1])) & 
                  (as.numeric(phecode) <= as.numeric(this_excluded_range_bounds[2])))
        
        #Remember that each exclusion phecode has a list of ICD10 codes that we should ignore.
        #Get the ICD10 or ICD9 codes which correspond to all the exclusion phecodes.
        this_exclude_icds = unlist(excluded_phecodes[[include_icd_column_name]])
        
        #Merge this exclusion criteria' ICD10 or ICD9 codes to all exclusion ICD10 or ICD9 codes for this phecode.
        exclude_icds = c(exclude_icds, this_exclude_icds) 
        
    }
    
    #Collapse the exclude icd10 or ICD9 list into a string to make it compatible with apply function.
    exclude_icds = paste(exclude_icds, collapse = ",")
        
    return(exclude_icds)
}

#Get the exclusion criteria for ICD10
phecode_map_aggregated = phecode_map_aggregated %>% 
    rowwise() %>% 
    mutate(ExcludeICD10s = getExclusionICD10s(phecode, IncludeICD10s, ExcludePheCodeRange, "IncludeICD10s"))

#The new column is a String, let's convert that into a list.
phecode_map_aggregated$ExcludeICD10s <- lapply(strsplit(
    as.character(phecode_map_aggregated$ExcludeICD10s), ","), trimws)

#Get the exclusion criteria for ICD9
phecode_map_aggregated = phecode_map_aggregated %>% 
    rowwise() %>% 
    mutate(ExcludeICD9s = getExclusionICD10s(phecode, IncludeICD9s, ExcludePheCodeRange, "IncludeICD9s"))

#The new column is a String, let's convert that into a list.
phecode_map_aggregated$ExcludeICD9s <- lapply(strsplit(
    as.character(phecode_map_aggregated$ExcludeICD9s), ","), trimws)

head(phecode_map_aggregated,50)

#Collpse
phecode_map_aggregated_flattened = phecode_map_aggregated %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = '|'))

write.table("Data/PhecodeFiles/PhecodeToICD.tsv", x = phecode_map_aggregated_flattened, row.names = F, sep = "\t", quote = F)

#How to read the phecode file back with columns as a list?
df = as.data.frame(fread("Data/PhecodeFiles/PhecodeToICD.tsv"))
phecode_map_aggregated = df %>% mutate_if(~any(str_detect(., fixed('|'))), ~str_split(., fixed('|')))
head(phecode_map_aggregated, 10)

phecode_map_aggregated[phecode_map_aggregated$phecode == "250.1",]

#Let's get df with all sample names with ICD10cm info.
aou_samples_with_icds = unique(disease_cohort_with_icds$person_id)

length(aou_samples_with_icds)

output_directory = "/home/jupyter/workspaces/geneticancestry/Data/PhecodeFiles/NewCohorts/"

createPheCodeFile = function(phecode_map_row_index){
    #Extract info from the phecode map df.
    phecode = phecode_map_aggregated[phecode_map_row_index,]$phecode
    
    ###Get the inclusion and exclusion criteria for ICD10 & ICD9.
    #We have to use [[1]] to extract the list out of the df.
    include_icd10s = phecode_map_aggregated[phecode_map_row_index,]$IncludeICD10s[[1]]
    exclude_icd10s = phecode_map_aggregated[phecode_map_row_index,]$ExcludeICD10s[[1]]

    #We extract for ICD9 now.
    include_icd9s = phecode_map_aggregated[phecode_map_row_index,]$IncludeICD9s[[1]]
    exclude_icd9s = phecode_map_aggregated[phecode_map_row_index,]$ExcludeICD9s[[1]]

    
    #Create a placeholder cohort data frame.
    cohort = data.frame("PersonID" = aou_samples_with_icds, 
                        "Status" = rep("Control", length(aou_samples_with_icds)))    

    #Get the samples which are a case. We will call them confident cases.
    confident_cases = disease_cohort_with_icds %>% 
        select(person_id, condition_source_value, condition_start_datetime) %>%
        dplyr::filter((condition_source_value %in% include_icd10s) | 
                     (condition_source_value %in% include_icd9s)) %>%
        select(person_id) %>% distinct(person_id) %>% pull(person_id)

    #Mark the cases.
    cohort = cohort %>% mutate(Status = ifelse(PersonID %in% confident_cases, "Case", Status))
    
    #For dubious codes, we will get those exclusion ICD10/ICD9 codes which are not present in include ICD10/ICD9 codes.
    #These samples are excluded from the cohort as such that they are nither control nor case.
    #We need to get ICD10/ICD9 codes that are present in Exclusion criteria, but not in inclusion criteria.
    #We are not confident if these are proper controls are not. They need to be excluded from the cohort.
    
    #We need to find codes that are in exclusion list but not in inclusoion list.
    dubious_icd10_codes = setdiff(exclude_icd10s, include_icd10s)
    dubious_icd9_codes = setdiff(exclude_icd9s, include_icd9s)
    
    dubious_samples = disease_cohort_with_icds %>% 
        select(person_id, condition_source_value) %>%
        filter((condition_source_value %in% dubious_icd10_codes) |
              (condition_source_value %in% dubious_icd9_codes)) %>%
        select(person_id) %>% distinct(person_id) %>% pull(person_id)
       
    #Mark the cases.
    cohort = cohort %>% mutate(Status = ifelse((PersonID %in% dubious_samples & Status == "Control"), "Excluded", Status))
        
    print(paste("Writing PheCode file for PheCode: ", phecode, sep = " "))
    
    #Write to a file for case control cohorts.
    output_file_name = paste("cohort",phecode,"tsv",sep = ".")
    
    
    #Write it to files, uncomment for final changes.
    write.table(x = cohort, file = paste(output_directory, output_file_name, sep = "/"),
                row.names = F, sep = "\t", quote = F)

}

#apply(phecode_map_aggregated[222:223,], 1, createPheCodeFile)

#Create phecode files now.
phecode_map_aggregated_indices = row.names(phecode_map_aggregated)
mclapply(phecode_map_aggregated_indices, function(i) createPheCodeFile(i), mc.cores = 6)

output_directory = "/home/jupyter/workspaces/geneticancestry/Data/PhecodeFiles/DoubleRuleCohorts//"


createPheCodeFileDoubleRule = function(phecode_map_row_index){
    #Extract info from the phecode map df.
    phecode = phecode_map_aggregated[phecode_map_row_index,]$phecode
    
    print(paste("Starting with Phecode:", phecode))
    ###Get the inclusion and exclusion criteria for ICD10 & ICD9.
    #We have to use [[1]] to extract the list out of the df.
    include_icd10s = phecode_map_aggregated[phecode_map_row_index,]$IncludeICD10s[[1]]
    exclude_icd10s = phecode_map_aggregated[phecode_map_row_index,]$ExcludeICD10s[[1]]

    #We extract for ICD9 now.
    include_icd9s = phecode_map_aggregated[phecode_map_row_index,]$IncludeICD9s[[1]]
    exclude_icd9s = phecode_map_aggregated[phecode_map_row_index,]$ExcludeICD9s[[1]]

    
    #Create a placeholder cohort data frame.
    cohort = data.frame("PersonID" = aou_samples_with_icds, 
                        "Status" = rep("Control", length(aou_samples_with_icds)))    

    #Get the samples which are a case. Count the number of assesments each sample has for a given date.
    cases = disease_cohort_with_icds %>% 
        select(person_id, condition_source_value, condition_start_datetime, DateOfAssesment) %>%
        dplyr::filter((condition_source_value %in% include_icd10s) | 
                     (condition_source_value %in% include_icd9s)) %>%
        #filter(person_id == "1043758" | person_id == "3488755") %>%    
        group_by(person_id) %>%
        count(DateOfAssesment, name= "NumberOfAssesments") 
        
    cases_double_rule = cases %>% 
                        count(person_id, name = "UniqueAssesments") %>% 
                        filter(UniqueAssesments > 1) %>%
                        select(person_id) %>% distinct(person_id) %>% 
                        pull(person_id)

    print(paste("Number of cases", length(cases_double_rule)))
    
    #Mark the cases.
    cohort = cohort %>% mutate(Status = ifelse(PersonID %in% cases_double_rule, "Case", Status))
    
    
    #For dubious codes, we will get those exclusion ICD10/ICD9 codes which are not present in include ICD10/ICD9 codes.
    #These samples are excluded from the cohort (only controls, not cases) 
    #as such that they are not a healthy control/
    #We need to get ICD10/ICD9 codes that are present in Exclusion criteria, but not in inclusion criteria.
    #We are not confident if these are proper controls are not. They need to be excluded from the cohort.
    
    #We need to find codes that are in exclusion list but not in inclusoion list.
    dubious_icd10_codes = setdiff(exclude_icd10s, include_icd10s)
    dubious_icd9_codes = setdiff(exclude_icd9s, include_icd9s)
    
    dubious_samples = disease_cohort_with_icds %>% 
        select(person_id, condition_source_value) %>%
        filter((condition_source_value %in% dubious_icd10_codes) |
              (condition_source_value %in% dubious_icd9_codes)) %>%
        select(person_id) %>% distinct(person_id) %>% pull(person_id)
       
    #Mark the cases.
    cohort = cohort %>% mutate(Status = ifelse((PersonID %in% dubious_samples & Status == "Control"), "Excluded", Status))
        
    print(paste("Writing PheCode file for PheCode: ", phecode, sep = " "))
    
    #Write to a file for case control cohorts.
    output_file_name = paste("cohort",phecode,"tsv",sep = ".")
    write.table(x = cohort, file = paste(output_directory, output_file_name, sep = "/"),
                row.names = F, sep = "\t", quote = F)


}

#apply(phecode_map_aggregated[222:223,], 1, createPheCodeFile)

#Create phecode files now.
phecode_map_aggregated_indices = row.names(phecode_map_aggregated)
mclapply(phecode_map_aggregated_indices, function(i) createPheCodeFileDoubleRule(i), mc.cores = 6)



criteria = c("E10.630", "E10.618", "E10.62" , "E10.9" , "E10.628", "E10.65" , "E10.621", "E10.8" , "E10.69" , 
             "E10.649", "E10.610", "E10.620", "E10" , "E10.622", "E10.6" , "E10.638", "E10.641", 
             "E10.61" , "E10.64" , "E10.63", "250.01", "250.03", "250.11", "250.13", "250.21", "250.23", 
             "250.31", "250.33", "250.41", "250.43", "250.51", "250.53", "250.61", "250.63", "250.71", "250.73", 
             "250.81", "250.83", "250.91", "250.93")
temp = disease_cohort_with_icds %>% filter(source_concept_code %in% criteria)




