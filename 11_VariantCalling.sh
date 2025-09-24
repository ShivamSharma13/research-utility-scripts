#!/bin/bash

# Place your gatech userame in the below export
export NAME="ssharma454"

get_input () {
	# Function for doing your getopts
	# Input: Getopts array
	while getopts "a:b:r:f:o:ezvih" option;
	do
        	case ${option} in
                	a) reads1=$OPTARG;;
                	b) reads2=$OPTARG;;
                	r) ref=$OPTARG;;
                	f) millsFile=$OPTARG;;
                	o) output=$OPTARG;;
			e) realign=1;;
                	z) gunzip=1;;
                	v) v=1;;
                	i) index=1;;
                	h) info_usage=1;;
			*) echo "Option not found. Please try again."
        	esac
	done

	if (( info_usage  ))
		then
			echo "a: path to reads file 1, b: path to reads file 2, r: reference file, f: millsFile for indels [use b38 version] , o: output file path, e: if you want to realign, z: if you want to zip the output file, i: if you want to index the imporvement file."
			exit 0 
	fi
	
	if (( v ))
	then
		:
	fi

	return 0
}

check_files () {
	# Function for checking for presence of input files, reference genome,
	# and the output VCF file
	#
	# Input: File locations (string)
	# Output: True, if checks pass; False, if checks fail (bool)
	
	#Checking
	
	#See if input was given or not.

	if (( v ))
	then
		echo "Checking for essential files for pipeline to run."
	fi

	if [ "$ref" == "" ] || [ "$reads1" == "" ] || [ "$reads2" == "" ] || [ "$millsFile" == "" ]
		then
			echo "Please send the complete arguments. "
			return 1
	fi
	
	#Checking for chromosome reference.
	file_check=0
	if [ ! -e "$ref" ];
		then
			echo "Chromosome reference RefSeq file not found."
			file_check=1
	fi
	
	#Checking for read file 1.
	if [ ! -e "$reads1" ];
       		then
                	echo "D2-DS3_paired1.fq file not found."   
			file_check=1
        fi
	
	#Checking for read file 2.
	if [ ! -e "$reads2" ];
        	then
                	echo "D2-DS3_paired2.fq file not found." 
			file_check=1
        fi

	#Checking for mills file.
	if [ ! -e "$millsFile" ];
		then
			echo "millsFile not found."
			file_check=1
	fi


	#Checking if output files exist.
	if [ -e "$output" ];
                then
			echo "Output file is already present. Do you want to overwite the existing file or exit? (y/n):"
			read -r choice
			if [ "$choice" == "y" ]
				then
					:
			elif [ "$choice" == "n" ]
				then	
					exit 1
			fi
        fi
	return $file_check
}

prepare_temp () {

	if (( v ))
	then
		echo "Setting up path variables."
	fi
	# Preparing your temporary directory
	current_path=$(pwd)
	tmp_path="${current_path}/tmp"
	data_path="${current_path}/data"
	lib_path="${current_path}/lib"
	if [ -e "${ref}.bwt" ]
	then
		echo "Found index file, skipping indexing step."
		return 0
	fi

	if (( v ))
	then
		echo "Creating index file from reference file using BWA."
	fi
	
	bwa index "$ref"	

	return 0
}


mapping () {
	# Function for the mapping step of the SNP-calling pipeline
	#
	# Input: File locations (string), Verbose flag (bool)
	# Output: File locations (string)
		
	#Creating lane.bam file.

	if (( v ))
	then
		echo "Creating lanes file using reference and read file."
	fi

	if [ -e "${tmp_path}/lane.sam" ]
	then
		echo "lane.sam found in tmp directory, using that..."
	else
		bwa mem -R '@RG\tID:chr17id\tSM:tmp\tLB:lib' "$ref" "$reads1" "$reads2" > "${tmp_path}/lane.sam"
	fi
	

	if (( v ))
	then
		echo "Succesfully created the lanes file, will run fixmate and sort on the result now..."
	fi

	#Imporving the lane.sam file
	if [ -e "${tmp_path}/lane_fixmate.bam" ]
	then
		echo "fixmate file found as well, using that..."
	else			
		samtools fixmate -O bam "${tmp_path}/lane.sam" "${tmp_path}/lane_fixmate.bam"
	fi
	

	#Sorting
	if [ -e "${tmp_path}/lane_sorted.bam" ]
	then
		echo "sorted lane file found, using that..."
	else
		samtools sort -O bam -o "${tmp_path}/lane_sorted.bam" -T "${tmp_path}/lane_temp" "${tmp_path}/lane_fixmate.bam"
	fi


	if (( v ))
	then
		echo "Successfuly created fixmate and sorted bam files."
	fi

	return 0
}

improvement () {
	# Function for improving the number of miscalls
	#
	# Input: File locations (string)
	# Output: File locations (string)
	
	#Preparing files required for imporving alignment.


	if (( v ))
	then
		echo "Running improvement functions now..."
	fi

	if [ -e "${tmp_path}/chr17.fa.fai" ]
		then
			echo "found an fa.fai file, using it..."
	else
		samtools faidx "$ref"
		samtools dict "$ref" -o "${data_path}/chr17.dict"
		
		samtools index "${tmp_path}/lane_sorted.bam"
	fi

	if (( v ))
	then
		echo "Success! Will run GATK aligners if asked for..."
	fi

	if (( realign  ))
	then	
		#Creating intervals file.
		if [ -e "${tmp_path}/lane.intervals" ]
			then
				echo "found lane.intervals file, using it..."
		else
			java -Xmx2g -jar "${lib_path}/GenomeAnalysisTK.jar" -T RealignerTargetCreator -R "$ref" -I "${tmp_path}/lane_sorted.bam" -o "${tmp_path}/lane.intervals" --known "$millsFile"		
		fi

		#Creating realigned file.
		if [ -e "${tmp_path}/lane_realigned.bam" ]
			then
				echo "found lane_realigned.bam, using it..."
		else
			java -Xmx4g -jar "${lib_path}/GenomeAnalysisTK.jar" -T IndelRealigner -R "$ref" -I "${tmp_path}/lane_sorted.bam" -targetIntervals "${tmp_path}/lane.intervals" -known "$millsFile" -o "${tmp_path}/lane_realigned.bam"
		
		fi


		if (( v ))
		then
			echo "Success! Will call variants now..."
		fi
	fi

	#Using base recalibrator
	#if [ -e "${tmp_path}/lane_recal.table" ]
	#then
		#echo "found lane_recal.table, using it..."
	#else
		:
		#java -Xmx4g -jar "${lib_path}/GenomeAnalysisTK.jar" -T BaseRecalibrator -R "$ref" -knownSites >bundle/b38/dbsnp_142.b38.vcf> -I <lane.bam> -o <lane_recal.table>
	#fi
	
	if (( index ))
	then
		samtools index "${tmp_path}/lane_realigned.bam"
	fi
	
	return 0
}

call_variants () {
	
	# Function to call variants
	#
	# Input: File locations (string)
	# Ouput: None
	
	if (( realign ))
	then
		align_file="${tmp_path}/lane_realigned.bam"
	else
		align_file="${tmp_path}/lane_sorted.bam"
	fi

	
	if (( v ))
	then
		echo "Will use bcftools for piling up files and call variants between them..."
	fi

	if (( gunzip ))
	then
			bcftools mpileup -Ob -o "${tmp_path}/result.bcf" -f "$ref" "${align_file}"; 
			bcftools call -vmO z -o "${output}.vcf.gz" "${tmp_path}/result.bcf";	
	else
		
			bcftools mpileup -Ob -o "${tmp_path}/result.bcf" -f "$ref" "${align_file}";
			bcftools call -vmO v -o "${output}.vcf" "${tmp_path}/result.bcf";
	fi


	if (( v ))
	then
		echo "Sucessfully created an output file."
		echo "Pipeline complete!"
		echo "Bye"
	fi

	return 0	
}


main() {
	# Function that defines the order in which functions will be called
	# You will see this construct and convention in a lot of structured code.
	
	# Add flow control as you see appropriate
	get_input "$@"
	
	if ! check_files "$1";
		then
			echo "Exiting because file check failed. Look for errors above."
			exit 1
	fi
	prepare_temp
	mapping # Add arguments here
	improvement # Add arguments here
	call_variants # Add arguments here
}

# Calling the main function
main "$@"


# DO NOT EDIT THE BELOW FUNCTION
bats_test (){
    command -v bats
}
# DO NOT EDIT THE ABOVE FUNCTION
