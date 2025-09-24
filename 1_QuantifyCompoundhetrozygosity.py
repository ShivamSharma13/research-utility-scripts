print("====> Importing HAIL, this might take a while...")
import hail as hl
import argparse
import pandas as pd
import os


'''
Define global variables pertinent to a dataset.
These variable are:
	variant_types_of_interest = Varinat labels from VEP indicating their level of pathogenicity.
'''

#variant_types_of_interest = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification']

variant_types_of_interest = ['missense_variant', 'protein_altering_variant', 'inframe_insertion', 'inframe_deletion']
genome_build = "GRCh38"

def prepare_cohort(vcf_file_path, fam_file_path, filter, filter_by):
	'''
	This function will read in the Phased VCF file and then merge it with the FAM file. 
	Merging with FAM file will let us filter any samples which are not of any interest to us.
	'''
	global genome_build
	
	print("====> Reading in the VCF file now using HAIL")
	if not os.path.exists(vcf_file_path):
		print("====> VCF file does not exist.")
		return -9

	#Read in the VCF file.
	recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
	phased_vcf = hl.import_vcf(vcf_file_path, reference_genome = genome_build, force_bgz = True, contig_recoding = recode)

	#Describe the VCF file.
	print("====> Number of samples in the VCF file = {}".format(phased_vcf.aggregate_cols(hl.struct(Records=hl.agg.count()))))

	#Show sample names from the VCF file.
	print("====> Here is what the VCF file looks like:")
	phased_vcf.s.show(5)
	phased_vcf.GT.show(5)
	
	
	print("====> Reading in the FAM file now using HAIL")
	#Read in the FAM file to get annotate and select only Cauc.
	fam_file = hl.import_table(fam_file_path, impute=True).key_by('Sample')
	fam_file.show(5)
	
	print("====> Annotating the VCF file using the FAM file now.")
	#Annotate the VCF object with the FAM file values.
	phased_vcf = phased_vcf.annotate_cols(pheno = fam_file[phased_vcf.s])
	phased_vcf.col.pheno.describe()
	
	#Look at the counts by Pop variable.
	print("====> Frequency of samples with disease status & Pop: ")
	print(dict(phased_vcf.aggregate_cols(hl.struct(pop_counts=hl.agg.counter(phased_vcf.col.pheno.Disease)))))
	print(dict(phased_vcf.aggregate_cols(hl.struct(pop_counts=hl.agg.counter(phased_vcf.col.pheno.Pop)))))

	if filter:
		#Filter out samples where Pop=Cauc
		phased_vcf_filtered = phased_vcf.filter_cols(phased_vcf.col.pheno.Pop == filter_by)

		print("====> Total number of Selected samples = {}".format(
			phased_vcf_filtered.aggregate_cols(hl.struct(pop_counts=hl.agg.counter(phased_vcf_filtered.col.pheno.Pop)))
			)
		)

		return phased_vcf_filtered
	
	else:
		return phased_vcf	  
	
def prepare_annotation_file(annotation_file_path, mane_file_path, genes_of_interest):
	'''
	This function will read in the VEP annotatio file and filter out non-connonical transcripts.
	Then it will select variants based on their pathogenicity.

	'''
	global variant_types_of_interest

	print("====> Reading in the VEP file now using HAIL")

	#Read the VEP file.
	vep_annotations = hl.import_table(annotation_file_path, comment = "^##.+", key = 'Feature').rename({'#Uploaded_variation': 'Uploaded_variation'})
	print("====> VEP annotation file looks like: ")
	print(vep_annotations.show(5))

	print("====> Reading in the MANE file now using HAIL")
	#Read the MANE cannonical transcripts file.
	cannonical_transcripts_mane = hl.import_table(mane_file_path, key = 'Ensembl_nuc').rename({'#NCBI_GeneID': 'NCBI_GeneID'})

	#Filter by genes of interest.
	cannonical_transcripts_mane = cannonical_transcripts_mane.filter(hl.set(genes_of_interest).contains(cannonical_transcripts_mane.symbol))

	#Show the cannonical_transcripts_mane filtered file now.
	print("====> MANE cannonical transcripts file looks like: ")
	print(cannonical_transcripts_mane.show(5))

	#Add a new column with version numbers of ENST IDs removed. Version numbers are absent in VEP annotations.                                                                   
	cannonical_transcripts_mane = cannonical_transcripts_mane.annotate(EnsembleIDNoVer = cannonical_transcripts_mane.Ensembl_nuc.replace('\.\d+', ''))

	#Re-key the table
	cannonical_transcripts_mane = cannonical_transcripts_mane.key_by('EnsembleIDNoVer')

	'''
	Get only the annotation rows of interest from VEP annotation.
	VEP annotation woul often times annotate variants with multiple transcript IDs
	We have decided to use only the connonical transcipts for this analyis.
	'''
	vep_annotations_select_transcript = vep_annotations.semi_join(cannonical_transcripts_mane)

	#Summarize the selected VEP annotations.
	vep_annotations_select_transcript_summary = vep_annotations_select_transcript.group_by(
		vep_annotations_select_transcript.Feature, vep_annotations_select_transcript.SYMBOL).aggregate(
		NumberOfVariants = hl.agg.count()
	)
	print("====> Variants present in each selected transcript: ")
	print(vep_annotations_select_transcript_summary.show())


	#Check the counts of different types of variants.
	#VEP gives a frozen_dict for aggregate function, we will overlap it with the LIST defined 
	#above to get all possible annotations.
	#Case example: The consequence "missense_variant,splice_region_variant" 
	#should also be included when looking for missense_variant.

	vep_variant_annotations = list(dict(dict(vep_annotations_select_transcript.aggregate(
		hl.struct(ConsequenceCount = hl.agg.counter(
			vep_annotations_select_transcript.Consequence)
				 )
	))['ConsequenceCount']).keys())

	print("====> Frequency table for variants: {}".format(vep_variant_annotations))

	#Get VEP annotation entries with at least one lof entries.
	vep_variant_annotations_of_interest = [variant_annotation for variant_type in variant_types_of_interest for variant_annotation in vep_variant_annotations if variant_type in variant_annotation]

	print("====> We will be looking at the following variants from the VEP file: {}".format(vep_variant_annotations_of_interest))

	#Get the loss of function variants. Use the list defined: 'vep_variant_annotations_of_interest' to filter these variants.
	vep_annotations_select_transcript_lof = vep_annotations_select_transcript.filter(
		hl.set(vep_variant_annotations_of_interest).contains(vep_annotations_select_transcript.Consequence))

	#Create a locus column similar to the MatrixTable.
	vep_annotations_select_transcript_lof = vep_annotations_select_transcript_lof.annotate(
		locus = "chr" + vep_annotations_select_transcript_lof.Location + ":" +
		vep_annotations_select_transcript_lof.REF_ALLELE + ":" +
		vep_annotations_select_transcript_lof.Allele
	)

	#Re-key the vep annotation df by locus.
	#We will use this annotation table to get only harmful variants from phased dataset.
	vep_annotations_select_transcript_lof = vep_annotations_select_transcript_lof.key_by(
		**hl.parse_variant(vep_annotations_select_transcript_lof.locus, reference_genome='GRCh38')
	)

	#Look at the table again.
	print("====> VEP variant table after Variant pthogenicity filtering now looks like: ")
	#Summarize the selected VEP annotations.
	vep_annotations_select_transcript_lof_summary = vep_annotations_select_transcript_lof.group_by(
		vep_annotations_select_transcript_lof.Feature, vep_annotations_select_transcript_lof.SYMBOL).aggregate(
		NumberOfVariants = hl.agg.count()
	)
	print(vep_annotations_select_transcript_lof_summary.show())


	return vep_annotations_select_transcript_lof


def process_phased_vcf_for_selected_variants(phased_vcf_filtered, vep_annotations_select_transcript_lof):
	#Select only the variants which are present in thr VEP anotation LoF file.
	print("====> Variants present originally: {}".format(phased_vcf_filtered.count()))
	phased_vcf_filtered_selected_variants = phased_vcf_filtered.semi_join_rows(vep_annotations_select_transcript_lof)

	print("====> Variants present after merging using VEP annotations: {}".format(phased_vcf_filtered_selected_variants.count()))

	#Annotate the variants with GeneName, GeneID, and TranscriptID
	phased_vcf_filtered_selected_variants = phased_vcf_filtered_selected_variants.annotate_rows(
		GeneSymbol = vep_annotations_select_transcript_lof[
			phased_vcf_filtered_selected_variants.locus, phased_vcf_filtered_selected_variants.alleles
		].SYMBOL,
		GeneID = vep_annotations_select_transcript_lof[
			phased_vcf_filtered_selected_variants.locus, phased_vcf_filtered_selected_variants.alleles
		].Gene,
		TranscriptID = vep_annotations_select_transcript_lof[
			phased_vcf_filtered_selected_variants.locus, phased_vcf_filtered_selected_variants.alleles
		].Feature
	)

	return phased_vcf_filtered_selected_variants

def calculate_compound_hetrozygosity(phased_vcf_filtered_selected_variants):
	'''
	Use a bunch of aggregations.
	'''
	print("====> Now calculating compund hetrozygosity numbers for each sample...")

	phased_vcf_aggregated_hetrozygosity_by_trascript = (phased_vcf_filtered_selected_variants.group_rows_by(
		phased_vcf_filtered_selected_variants.TranscriptID)
		.aggregate_entries(Hetrozygous1AND0 = hl.agg.count_where(
			(phased_vcf_filtered_selected_variants.GT[0] == 1) & (phased_vcf_filtered_selected_variants.GT[1] == 0)
		),
			Hetrozygous0AND1 = hl.agg.count_where(
			(phased_vcf_filtered_selected_variants.GT[0] == 0) & (phased_vcf_filtered_selected_variants.GT[1] == 1)
		))
		.result()
	)

	print("====> Variants present after merging using VEP annotations: ")
	print(phased_vcf_aggregated_hetrozygosity_by_trascript.show(n_cols = 20))

	phased_vcf_aggregated_compound_hetrozygous = (phased_vcf_aggregated_hetrozygosity_by_trascript.group_rows_by(
	phased_vcf_aggregated_hetrozygosity_by_trascript.TranscriptID)
		.aggregate(IsCompoundHetrozygous = hl.agg.count_where(
			(phased_vcf_aggregated_hetrozygosity_by_trascript.Hetrozygous0AND1 > 0) 
			& 
			(phased_vcf_aggregated_hetrozygosity_by_trascript.Hetrozygous1AND0 > 0))
		)
	)

	print("====> Hetrozygosity counts for each sample by transcript ID: ")
	print(phased_vcf_aggregated_compound_hetrozygous.show(n_cols = 20))

	return phased_vcf_aggregated_compound_hetrozygous


def dump_compound_hetrozygous_numbers(phased_vcf_aggregated_compound_hetrozygous, output_file_path):
	print("====> Now dumping results...")

	#Export the aggregated file to a temporary TSV file.
	phased_vcf_aggregated_compound_hetrozygous.IsCompoundHetrozygous.export(".tmp.IsCompoundHetrozygous.tsv")

	#Read the file back in.
	current_aggregated_compound_hetrozygous_numbers = pd.read_csv(".tmp.IsCompoundHetrozygous.tsv",sep="\t")

	current_aggregated_compound_hetrozygous_numbers = current_aggregated_compound_hetrozygous_numbers.transpose()

	#Adjust index and col names.
	current_aggregated_compound_hetrozygous_numbers = current_aggregated_compound_hetrozygous_numbers.reset_index(level=0) # Set index to a column
	current_aggregated_compound_hetrozygous_numbers = current_aggregated_compound_hetrozygous_numbers.rename(
		columns=current_aggregated_compound_hetrozygous_numbers.iloc[0]
		).drop(current_aggregated_compound_hetrozygous_numbers.index[0]) # Move 1st row to header and drop row with header

	#Check if output file already exists.
	if os.path.exists(output_file_path):
		#Merge the columns together.		
		existing_aggregated_compound_hetrozygous_numbers = pd.read_csv(output_file_path, sep="\t", header=0, index_col=0)
		combined_aggregated_compound_hetrozygous_numbers = pd.merge(existing_aggregated_compound_hetrozygous_numbers,current_aggregated_compound_hetrozygous_numbers,on="TranscriptID",how="outer") # Outer join to place NA but maintain all samples

		#print(combined_aggregated_compound_hetrozygous_numbers)
		combined_aggregated_compound_hetrozygous_numbers.to_csv(output_file_path, sep = "\t")
		
	else:
		#Create this file for the first time.
		current_aggregated_compound_hetrozygous_numbers.to_csv(output_file_path, sep = "\t")

	return

def main():
	'''
	Main function will be reading:
		1. FAM File indicating the sample names, and if they need to be filtered.
		2. VEP annotation file.
		3. Phased VCF file for the cohort.
		4. MANE file which has the information for the cannonical trancripts.
		5. Genes of interest file.
	'''

	parser = argparse.ArgumentParser()

	#Arguments.
	parser.add_argument("-f", "--fam", help="Path to FAM file. FAM file should have: 1 = Sample, 6 = Disease, 7 = SamplePop", required=True)
	parser.add_argument("-a", "--annotation", help="Path to VEP Annotation file.", required=True)
	parser.add_argument("-m", "--mane", help="Path to the MANE file.", required = True)
	parser.add_argument("-v", "--vcf", help="Path to phased VCF file", required=True)
	parser.add_argument("-g", "--gene", help="Path to genes of interest file", required=True)
	parser.add_argument("-o", "--output", help="Path to output file", required=True)


	
	#Parse and gather whatever the user sent.
	args = vars(parser.parse_args())
	fam_file_path = args['fam']
	annotation_file_path = args['annotation']
	mane_file_path = args['mane']
	vcf_file_path = args['vcf']
	genes_of_interest_file_path = args['gene']
	output_file_path = args['output']
	

	#Obtain the genes of interest. This file will have one gene symbol per line.
	with open(genes_of_interest_file_path) as f:
		raw = f.read()

	genes_of_interest = [i for i in raw.split("\n") if i != ""]
	
	print("====> Genes of interest file read. Here are the {} genes read: {}".format(len(genes_of_interest), genes_of_interest))
	
	#Initialize HAIL.
	print("====> Now initializing HAIL, this might take a while...")
	hl.init()
	

	print("====> Now preparing cohort...")
	#Prepare the cohort using HAIL.
	phased_vcf_filtered = prepare_cohort(vcf_file_path, fam_file_path, filter = True, filter_by = 'Cauc')

	if phased_vcf_filtered == -9:
		print("====> Skipping this chromosome...")
		return

	vep_annotations_select_transcript_lof = prepare_annotation_file(annotation_file_path, mane_file_path, genes_of_interest)

	phased_vcf_filtered_selected_variants = process_phased_vcf_for_selected_variants(phased_vcf_filtered, vep_annotations_select_transcript_lof)

	phased_vcf_aggregated_compound_hetrozygous = calculate_compound_hetrozygosity(phased_vcf_filtered_selected_variants)

	dump_compound_hetrozygous_numbers(phased_vcf_aggregated_compound_hetrozygous, output_file_path)


if __name__ == "__main__":
	main()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	