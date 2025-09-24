rm(list=ls())
set.seed(13)

options(scipen=100, digits=3)

library('dplyr')
library('ggplot2')

#Set the working directory.
setwd("Data/")

#Load the Demographic file and check it out.
demographics <- read.table('Metadata/DemographicData.tsv', header = T, sep = '\t')

#Get relevant columns.
demographics = demographics %>% select(person_id, race, ethnicity)

head(demographics)

#Read rye file.
rye_estimates <- read.table('Rye/extractedChrAllPruned.30.Rye-30.13.Q', header = T, sep = '\t')

#check the matrix.
head(rye_estimates)

#Merge the dataset.
ancestry_dataset = merge(demographics, rye_estimates, by.x = "person_id", by.y = "SampleID")

#Get some aggregate information.
print(paste("Total samples = ", nrow(ancestry_dataset), sep = " "))

print("Sample count breakdown:")
count(ancestry_dataset, race)

#Create a new variable with original race values except for Hisapanic/Latino folks, 
#who will get their race labeled as Hispanic/Latino.

#########################LOGIC#############################

###First, all the people who marked one of the none options in their race option will get "No information" labeled.
#Create a new variable with original race values except for Hisapanic/Lation folks, 
#who will get their race labeled as Hispanic/Latino.

ancestry_dataset$SelfReportedRaceEthnicity = ifelse((ancestry_dataset$race == "None Indicated" | 
                                                ancestry_dataset$race == "I prefer not to answer" | 
                                                ancestry_dataset$race == "None of these" | 
                                                ancestry_dataset$race == "PMI: Skip"), 
                                               "No information", ancestry_dataset$race)

count(ancestry_dataset, SelfReportedRaceEthnicity)


###Second, people who marked "Hispanic or latino" in ethnicity will get that label in SelfReportedRaceEthnicity.
ancestry_dataset$SelfReportedRaceEthnicity = ifelse((ancestry_dataset$ethnicity == "Hispanic or Latino" &
                                                ancestry_dataset$SelfReportedRaceEthnicity == "No information"), 
                                           "Hispanic or Latino", 
                                           ancestry_dataset$SelfReportedRaceEthnicity)   

count(ancestry_dataset, SelfReportedRaceEthnicity)


###Third, people who marked Hispanic or latino in their ethnicity but did not have "No information"
#in their race, will get "More than one population" label.
ancestry_dataset$SelfReportedRaceEthnicity = ifelse((ancestry_dataset$ethnicity == "Hispanic or Latino" &
                                                ancestry_dataset$SelfReportedRaceEthnicity != "Hispanic or Latino"), 
                                           "More than one population", 
                                           ancestry_dataset$SelfReportedRaceEthnicity)   

count(ancestry_dataset, SelfReportedRaceEthnicity)


#Multiply all estimates by 100.
multiply <- function(x, na.rm = FALSE) (x*100)

#Mutatet columns.
ancestry_dataset = ancestry_dataset %>% mutate_at(vars(SouthAsianIndian:NativeAmerican), multiply)

#Let's look at the dataset briefly.
head(ancestry_dataset)

#Collapse the 14 sub-gropus into 7 continental groups.
admixture_dataset_continental = ancestry_dataset 
admixture_dataset_continental$African = admixture_dataset_continental$AfricanNigerian +
                                admixture_dataset_continental$AfricanSeneGambian +
                                admixture_dataset_continental$AfricanBantu

admixture_dataset_continental$European = admixture_dataset_continental$EuropeanFinnish +
                                admixture_dataset_continental$EuropeanBritish +
                                admixture_dataset_continental$EuropeanSpanish +
                                admixture_dataset_continental$EuropeanItalian

admixture_dataset_continental$Asian = admixture_dataset_continental$SouthAsianPakistan + 
                                admixture_dataset_continental$SouthAsianIndian + 
                                admixture_dataset_continental$EastAsian

admixture_dataset_continental = admixture_dataset_continental %>% 
                                                        select(person_id, SelfReportedRaceEthnicity,
                                                       African, European, Asian,
                                                       NativeAmerican, WestAsian, Oceanian)

head(admixture_dataset_continental)

write.table(x = admixture_dataset_continental, file = "Rye/extractedChrAllPruned.ContinentalEstimates.6.Q", 
            row.names = F, sep = "\t", quote = F)

#Define colors for admixture plots.
continental_ancestry_colors = c("Asian" = "#B71C1C", "WestAsian" = "#80461B", 
                               "European" = "#FFA000", "African" = "#283593",
                               "NativeAmerican" = "#41C9F8", "Oceanian" = "#E040FB")

plot_admixture <- function(dataset){
    #Get the race for this dataset.
    this_dataset_race = unique(dataset$SelfReportedRaceEthnicity)

    #Get the mean ancestry value for each column.
    means = as.data.frame.list(colMeans(dataset %>% 
                                select(African, European, Asian, 
                                       NativeAmerican, WestAsian, Oceanian)))
    
    #Get the max ancestry column index.
    max_ancestry_column_index = max.col(means)
    
    #Get max ancestry column name.
    max_ancestry_column_name = colnames(means)[max_ancestry_column_index]
    
    #Sort the dataset by max ancestry column.
    dataset_sorted = dataset[order(dataset[,max_ancestry_column_name], decreasing = TRUE),]

    row.names(dataset_sorted) <- NULL
    
    #Melt the dataset and preserve the index for the admixture plot.
    dataset_sorted$index= as.numeric(rownames(dataset_sorted))

    dataset_sorted_melt = reshape2::melt(data = dataset_sorted, 
                                        id.vars = c('person_id', 'SelfReportedRaceEthnicity', 'index'), 
                                        measure.vars = c('European', 'African', 'NativeAmerican',
                                                   'Asian', 'Oceanian', 'WestAsian'))
    colnames(dataset_sorted_melt)[4] <- 'Ancestry'
    
    #Make the ggplot.
    rye_continental_admixture_plot = ggplot(data = dataset_sorted_melt,
                            aes(x=as.factor(index), y=value, fill=Ancestry)) +
    geom_bar(stat="identity", width=1) + 
    scale_fill_manual(values=continental_ancestry_colors, guide = 'none') + 
    labs(title = "") +
    theme(line = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), 
          axis.text=element_blank(), axis.ticks=element_blank(), panel.spacing.x=unit(1, "lines"), 
          strip.background = element_blank(), panel.background = element_blank(), 
          strip.text.x = element_text(size = 24), title = element_text(size = 32), 
          legend.text=element_text(size=24), legend.title = element_text(size = 28))

    #Save the plot.
    ggsave(paste("../Figures/Admixture", gsub(" ", "_", this_dataset_race, fixed = TRUE), "png", sep = "."),
           rye_continental_admixture_plot, width = 24, height = 16, dpi = 300, 
           units = "in", device='png', limitsize = FALSE)
    ggsave(paste("../Figures/Admixture", gsub(" ", "_", this_dataset_race, fixed = TRUE), "pdf", sep = "."),
           rye_continental_admixture_plot, width = 24, height = 16, dpi = 300, 
           units = "in", device='pdf', limitsize = FALSE)
    
    return(means)
}



#Remove the missing race individuals.
admixture_dataset_continental %>% 
   group_by(SelfReportedRaceEthnicity) %>%
   do(plot_admixture(.))


