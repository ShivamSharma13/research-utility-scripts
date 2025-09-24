# ResearchScripts
The repository has multiple scripts in R &amp; Python that I have used during my PhD to perform certain analysis in the broad area of human population genetics, health disparities, and NGS data analysis.

## Here is a brief summary of these scripts:

#### 1.  1_QuantifyCompoundhetrozygosity.py
This python based script uses HAIL to quantify events of compound hetrozygosity (where someone inherits two copies of a deleterious variants in a gene) in multi-sample genetic data. We can customize different masks (such as missense, frameshift, or stop gained) to filter for these variants. 

This script has been tested on WGS data of size ~3,000 samples and WES data of size ~50,000 samples. It finishes within an hour - thanks to HAIL!


#### 2. 2_ModelDiseaseLandscape.R
This R script performs a phenome wide regression testing of three ancestry estimates (European, African, & Native American) in a biobank-scale dataset. The whole phenome set of size ~1,700 is defined by Phecode classification which are harmonized using ICD10cm and ICD9cm codes from electronic health data of individuals in the biobank. 

This script is tested to work on biobanks with genetic and EHR data with sample sizes > 100,000


#### 3. 3_BedtoolsOverlap.py
A Python-based implementation of Bedtools' overlap functionality which will take fixed (and low!) memory. This aproach uses a sliding-window overlap matrix to find overlapping bedding file coordinates. This script is faster than the original tool.


#### 4. 4_RNAseqDiffrentialGeneExpression.py
A python pipeline for quantification of RNA-seq samples. This is the key step for an RNA-seq pipeline. I used [Salmon](https://combine-lab.github.io/salmon/) for this analysis.


#### 5. 5_MediationAnalysis.R
This script performs a mediation analysis across the whole genome to test wheather a copy of variant is mediating a phenotype (disease) of interest. The script works pretty similar to how a GWAS analysis works - it independently tests for each SNP to check whether it mediates (instead of association) a relationship between phenotype of interest (disease) and an associated phenotype (socio-economic deprivation).


#### 6. 6_CreateDosageMatrix.R
Sometimes it is benificial to get dosage values {0,1,2} from genetic data which can further be used downstream to analyze effects of certain variants on a phenotype of interest. Dosages denote how many copies of a variant a sample might have? 0, 1 or 2. This script uses data frame manipulation and plink commands to produce the output matrix.


#### 7. 7_CreatePhecodeCohorts.R

This script is useful for creating phecode cohorts to define samples as either Cases, Conrtrols, or Excluded. These are helpfule for performing association testing where cases and controls can be defined with high confidence. It can run on multiple cores with high speed.

#### 8. 8_All2FASTA.py
Too many formats for sequences but too little time! This script infers a sequencing filetype (fastq, embl, genbank or MEGA) and converts to FASTA format.

#### 9. 9_Transcript2GenomeMapper.py
NCBI coordinates for a matured mRNA NM_xxx like entry are nothing more than an index. On the other hand, popular miRNA binding site databases also use these coordinates instead of Genome Coordinates. For analysis of these binding sites for our research, we need genomic coordinates instead of the coordinates mentioned here. This script uses NCBI mRNAs and maps that mRNA's exon coordinates to the genome using the GENCODE annotation GFF3 file. 

#### 11. 11_VariantCalling.sh
A classic variant calling pipeline which uses bwa, samtools and gatk for simple variant calling.

#### 12. 12_PlotAdmixturePlots.R
Admixture plots are helpful to look at the genetic diversity of large-scale NGS/Genotyping datasets. Population structure is an essential confounder in association studies and needs to be accounted for in various statistical methods. This script generates admixture plots for diverse datasets stratified by the self-reported race/ethnicity information.



