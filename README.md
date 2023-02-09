# ResearchScripts
The repository has multiple scripts in R &amp; Python that I have used during my PhD to perform certain analysis in the broad area of human population genetics, health disparities and NGS data analysis.

## Here is a brief summary of these scripts:

#### 1.  quantifyCompoundhetrozygosity.py 
This python based script uses HAIL to quantify events of compound hetrozygosity (where someone inherits two copies of a deleterious variants in a gene) in multi-sample genetic data. We can customize different masks (such as missense, frameshift, or stop gained) to filter for these variants. 

This script has been tested on WGS data of size ~3,000 samples and WES data of size ~50,000 samples. It finishes within an hour - thanks to HAIL!

#### 2. modelDiseaseLandscape.R
This R script performs a phenome wide regression testing of three ancestry estimates (European, African, & Native American) in a biobank-scale dataset. The whole phenome set of size ~1,700 is defined by Phecode classification which are harmonized using ICD10cm and ICD9cm codes from electronic health data of individuals in the biobank. 

This script is tested to work on biobanks with genetic and EHR data with sample sizes > 100,000

#### 3. createPhecodeCohorts.R

This script is useful for creating phecode cohorts to define samples as either Cases, Conrtrols, or Excluded. These are helpfule for performing association testing where cases and controls can be defined with high confidence. It can run on multiple cores with high speed.

#### 4. plotAdmixturePlots.R

Admixture plots are helpful to look at the genetic diversity of large-scale NGS/Genotyping datasets. Population structure is an essential confounder in association studies and needs to be accounted for in various statistical methods. This script generates admixture plots for diverse datasets stratified by the self-reported race/ethnicity information.
