#!/usr/bin/env python3
import argparse

def line_to_values(line):
	#Function to split a line and send bed file values.
	line = line.rstrip('\n')
	chromosome_name = line.split('\t')[0]
	chromosome_start = int(line.split('\t')[1])
	chromosome_stop = int(line.split('\t')[2])
	misc_rest = line.split('\t')[3:]

	return chromosome_name, chromosome_start, chromosome_stop, misc_rest


def see_for_overlaps(input_file_1_path, input_file_2_path, minimal_overlap, join, output_file_path):
	file_1_chromosome_memory = None
	file_2_chromosome_memory = None

	file_2_handle_memory = []

	with open(input_file_1_path, 'r') as file_1_handle, open(input_file_2_path, 'r') as file_2_handle, open(output_file_path, 'w') as op_file_handle:
		file_1_chromosome_name, file_1_chromosome_start, file_1_chromosome_stop, file_1_misc_rest = line_to_values(file_1_handle.readline())
		if join:
			file_2_handle_memory.append(file_2_handle.tell())

		file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())
		

		#Set up memories.
		file_1_chromosome_memory = file_1_chromosome_name
		file_2_chromosome_memory = file_2_chromosome_name

		while True:
			try:
				#print(file_2_handle_memory)
				#print("Here" ,file_1_chromosome_start, file_1_chromosome_stop)

				#################################################################################
				########################Chromosome Change Logic##################################
				if file_1_chromosome_name != file_1_chromosome_memory and file_1_chromosome_memory == file_2_chromosome_name:
					#File 2 is behind.
					temp = file_2_handle.tell()
					file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())
					if file_2_chromosome_name == file_1_chromosome_name:
						file_1_chromosome_memory = file_1_chromosome_name
						file_2_chromosome_memory = file_2_chromosome_name
						file_2_handle_memory.clear()
						file_2_handle_memory.append(temp)
						continue
					else:
						continue

				if file_2_chromosome_name != file_2_chromosome_memory and file_2_chromosome_memory == file_1_chromosome_name:
					#File 1 is behind.
					file_1_chromosome_name, file_1_chromosome_start, file_1_chromosome_stop, file_1_misc_rest = line_to_values(file_1_handle.readline())
					if join:
						file_2_handle_memory.clear()
					if file_1_chromosome_name == file_2_chromosome_name:
						file_2_chromosome_memory = file_2_chromosome_name
						file_1_chromosome_memory = file_1_chromosome_name
						continue
					else:
						continue

				###########################Chromosome Change Logic ENDS##############################


				########################################################################################	
				##################One File's pointer is behind another file's pointer###################
				#########################File Pointer Shifting Logic Begins#############################

				if file_1_chromosome_name == file_2_chromosome_name and file_1_chromosome_stop <= file_2_chromosome_start:
					#Move the pointer of file 1.
					#print(file_1_chromosome_start, file_1_chromosome_stop)
					file_1_chromosome_name, file_1_chromosome_start, file_1_chromosome_stop, file_1_misc_rest = line_to_values(file_1_handle.readline())
					if join and len(file_2_handle_memory) > 0:
						#if Join is true then move the pointer back.
						file_2_handle.seek(file_2_handle_memory[0])
						file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())
						file_2_handle_memory.clear()
						#print(file_1_chromosome_start, file_1_chromosome_stop)

						#print(str([file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest]))
					continue

				if file_1_chromosome_name == file_2_chromosome_name and file_1_chromosome_start >= file_2_chromosome_stop:
					#Move the pointer of file 2.
					if join:
						file_2_handle_memory.append(file_2_handle.tell())

					file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())
					continue

				##########################File Pointer Shifting Logic Ends##############################


				if file_1_chromosome_name == file_2_chromosome_name and file_1_chromosome_start >= file_2_chromosome_start and file_1_chromosome_start <= file_2_chromosome_stop:
					if file_1_chromosome_stop <= file_2_chromosome_stop:
						############Complete Overlap.#################
						if join:
							output_line_1 = [str(file_1_chromosome_name), str(file_1_chromosome_start), str(file_1_chromosome_stop)] + [str(i) for i in file_1_misc_rest]
							output_line_2 = [str(file_2_chromosome_name), str(file_2_chromosome_start), str(file_2_chromosome_stop)] + [str(i) for i in file_2_misc_rest]
							
							op_file_handle.write("\t".join(output_line_1 + output_line_2) + '\n')
							file_2_handle_memory.append(file_2_handle.tell())
							file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())
							
							continue
						else:
							#print(str(file_1_chromosome_name) + "\t" + str(file_1_chromosome_start) + "\t" + str(file_1_chromosome_stop) + "\t" + "\t".join([str(i) for i in file_1_misc_rest]))
							output_line = [str(file_1_chromosome_name), str(file_1_chromosome_start), str(file_1_chromosome_stop)] + [str(i) for i in file_1_misc_rest]
							op_file_handle.write("\t".join(output_line) + '\n')
							file_1_chromosome_name, file_1_chromosome_start, file_1_chromosome_stop, file_1_misc_rest = line_to_values(file_1_handle.readline())
							continue	
						##########Complete Overlap Ends##############


					if file_1_chromosome_stop >= file_2_chromosome_stop:
						##############Partial Overlap###################
						file_1_entry_length = file_1_chromosome_stop - file_1_chromosome_start
						if (file_2_chromosome_stop - file_1_chromosome_start)/file_1_entry_length >= minimal_overlap/100:
							if join:
								output_line_1 = [str(file_1_chromosome_name), str(file_1_chromosome_start), str(file_1_chromosome_stop)] + [str(i) for i in file_1_misc_rest]
								output_line_2 = [str(file_2_chromosome_name), str(file_2_chromosome_start), str(file_2_chromosome_stop)] + [str(i) for i in file_2_misc_rest]
								
								op_file_handle.write("\t".join(output_line_1 + output_line_2) + '\n')
								file_2_handle_memory.append(file_2_handle.tell())
								file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())
								continue
							else:
								#print(str(file_1_chromosome_name) + "\t" + str(file_1_chromosome_start) + "\t" + str(file_1_chromosome_stop) + "\t" + "\t".join([str(i) for i in file_1_misc_rest]))
								output_line = [str(file_1_chromosome_name), str(file_1_chromosome_start), str(file_1_chromosome_stop)] + [str(i) for i in file_1_misc_rest]
								op_file_handle.write("\t".join(output_line) + '\n')
								file_1_chromosome_name, file_1_chromosome_start, file_1_chromosome_stop, file_1_misc_rest = line_to_values(file_1_handle.readline())
								continue
						else:
							#Move on to the next.
							file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())						
							continue
						############Partial Overlap Ends#################

				##############################Other Partial Overlap Logic#######################
				if file_1_chromosome_name == file_2_chromosome_name and file_1_chromosome_start <= file_2_chromosome_start and file_1_chromosome_stop >= file_2_chromosome_start and file_1_chromosome_stop <= file_2_chromosome_stop:
					file_1_entry_length = file_1_chromosome_stop - file_1_chromosome_start
					if (file_1_chromosome_stop - file_2_chromosome_start)/file_1_entry_length >= minimal_overlap/100:
						if join:
							output_line_1 = [str(file_1_chromosome_name), str(file_1_chromosome_start), str(file_1_chromosome_stop)] + [str(i) for i in file_1_misc_rest]
							output_line_2 = [str(file_2_chromosome_name), str(file_2_chromosome_start), str(file_2_chromosome_stop)] + [str(i) for i in file_2_misc_rest]
							
							op_file_handle.write("\t".join(output_line_1 + output_line_2) + '\n')
							file_2_handle_memory.append(file_2_handle.tell())
							
							file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())
							continue
						else:
							#print(str(file_1_chromosome_name) + "\t" + str(file_1_chromosome_start) + "\t" + str(file_1_chromosome_stop) + "\t" + "\t".join([str(i) for i in file_1_misc_rest]))
							output_line = [str(file_1_chromosome_name), str(file_1_chromosome_start), str(file_1_chromosome_stop)] + [str(i) for i in file_1_misc_rest]
							op_file_handle.write("\t".join(output_line) + '\n')
							file_1_chromosome_name, file_1_chromosome_start, file_1_chromosome_stop, file_1_misc_rest = line_to_values(file_1_handle.readline())
							continue
				
				##########################Other Partial Overlap Logic Ends#######################

				##########################Situation where File 1 read completely engulfs the File 2 read#########
				if file_1_chromosome_name == file_2_chromosome_name and file_2_chromosome_start >= file_1_chromosome_start and file_2_chromosome_stop < file_1_chromosome_stop:
					############But Check if the completely inside read still covers more than minimal overlap value########
					if (file_2_chromosome_stop - file_2_chromosome_start) >= (minimal_overlap/100) * (file_1_chromosome_stop - file_1_chromosome_start):
						if join:
							output_line_1 = [str(file_1_chromosome_name), str(file_1_chromosome_start), str(file_1_chromosome_stop)] + [str(i) for i in file_1_misc_rest]
							output_line_2 = [str(file_2_chromosome_name), str(file_2_chromosome_start), str(file_2_chromosome_stop)] + [str(i) for i in file_2_misc_rest]
							
							op_file_handle.write("\t".join(output_line_1 + output_line_2) + '\n')
							file_2_handle_memory.append(file_2_handle.tell())
							
							file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())
							continue
						else:
							output_line = [str(file_1_chromosome_name), str(file_1_chromosome_start), str(file_1_chromosome_stop)] + [str(i) for i in file_1_misc_rest]
							op_file_handle.write("\t".join(output_line) + '\n')
							file_1_chromosome_name, file_1_chromosome_start, file_1_chromosome_stop, file_1_misc_rest = line_to_values(file_1_handle.readline())
							continue
					##############Inner Segment Completely Inside but covers more than minimal overlap Ends##################
					file_2_chromosome_name, file_2_chromosome_start, file_2_chromosome_stop, file_2_misc_rest = line_to_values(file_2_handle.readline())
					continue

				file_1_chromosome_name, file_1_chromosome_start, file_1_chromosome_stop, file_1_misc_rest = line_to_values(file_1_handle.readline())
				continue
			except IndexError:
				break


def get_arguments():
	# Argparse code
	parser = argparse.ArgumentParser()
	parser.add_argument("-i1", "--input-file-1", help="Path to Input File 1.", required=True)
	parser.add_argument("-i2", "--input-file-2", help="Path to Input File 1.", required=True)
	parser.add_argument("-m", "--minimal-overlap", help="An Integer value for Minimal Overlap.", required=True)
	parser.add_argument("-j", "--join", help="Do you want to join the two entries?", required=False, action='store_true')
	parser.add_argument("-o", "--output-file", help="Path to Output File.", required=True)
	args = vars(parser.parse_args())

	input_file_1_path = args['input_file_1']
	input_file_2_path = args['input_file_2']
	minimal_overlap = int(args['minimal_overlap'])
	join = args['join']
	output_file_path = args['output_file']
	return input_file_1_path, input_file_2_path, minimal_overlap, join, output_file_path


if __name__ == "__main__":
	input_file_1_path, input_file_2_path, minimal_overlap, join, output_file_path = get_arguments()
	see_for_overlaps(input_file_1_path, input_file_2_path, minimal_overlap, join, output_file_path)




