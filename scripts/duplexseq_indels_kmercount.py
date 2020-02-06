#!/usr/bin/env python3

import sys
import csv
import re
import os
import copy

def indels_kmer():

	'''
	This function uses a kmer database(s) and a set of indels from our duplexseq processing pipeline to report kmer abundances
	
	It is called with:
	NameOfScript
		followed by these 5 arguments (all separated with spaces)
		1. input file with path
		2. kmer database path 
		3. list of kmer databases separated by commas (no spaces). For example: symb_shotgun_k39max10000_short.txt,symb_shotgun_k39max10000_short2.txt
		   if user only has one kmer database, that is fine, put spaces before and after and don't use any commas
		4. output file and path
		5. Total number of columns in input file
	
	
#example call: duplexseq_indels_kmercount.py  indels_file.txt  /home/KMC_databases/  my_kmer_db.txt  indels_file.kmer.txt  12
	'''

	inputfile = sys.argv[1]
	kmerdb_path = sys.argv[2]
	kmerdb_name_list_temp = sys.argv[3]
	outputfile = sys.argv[4]
	total_cols = int(sys.argv[5])

	comp_dict ={'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a'}
	
	###
	###
	###
	'''
	make a nested dictionary out of the kmer database(s)
	In the top level, each kmer db is a dictionary key. The value of each each key is another dictionary
	In this second, bottom level ditionary - each key is the kmer (the 39 nt sequence ) and each key value is the count
	'''
	###
	###
	###
	
	#open kmer files by making a list of all the files within a given replicate
	#the handle list is the same length as the files list but contains handle names instead of file names
	
	kmer_counter=0
	kmerdb_name_list=kmerdb_name_list_temp.split(',')
	kmer_handle_list = []
	kmer_topdict={}
	kmer_header=''
	
	for kmer_db_file in kmerdb_name_list:
		kmer_header = kmer_header + '\tAltcount_' + kmer_db_file
		kmer_counter += 1
		kmer_handle_list.append('file' + str(kmer_counter))
		kmer_topdict['kmer'+str(kmer_counter)] = {}
		try:
				kmer_handle_list[kmer_counter-1]=open(kmerdb_path + kmer_db_file)	
		except:
				return ('unable to open/find kmer db file(s): ' + str(kmerdb_path + kmer_db_file))
		
		kmer_temp_reader=csv.reader(kmer_handle_list[kmer_counter-1], delimiter = '\t')
		for row in kmer_temp_reader:
			kmer_topdict['kmer'+str(kmer_counter)][row[0]] =kmer_topdict['kmer'+str(kmer_counter)].get(row[0],row[1])
		
			
	indels_col_number = total_cols + len(kmerdb_name_list)
	
		
	###
	###
	###
	'''
	loop through the indels files, find kmers, and write to ouput
	'''
	###
	###
	###
	
	
	
	#create lists of file handles and output files for each user specified rep
	indels_counter=0
	
	#there is no longer a need for a loop here, but I am keeping the structure of gus's original code intact
	for dummy_var in [1]:
		
		indels_file_list=[]
		indels_handle_list=[]
		indels_reader_list=[]
		outer_low_list=[]
	
		indels_counter+=1
		indels_file_list.append(inputfile)
		indels_handle_list.append('file' + str(indels_counter))
		indels_reader_list.append(('reader' + str(indels_counter)))
		outer_low_list.append(outputfile)
		
		
		for h_num in range(len(indels_handle_list)):
			try:
				indels_handle_list[h_num]=open(indels_file_list[h_num])		
				outer_low_list[h_num]=open(outer_low_list[h_num], 'w')

			except:
				return ('unable to open/create the file: ' + indels_file_list[h_num] + ' or ' + outer_low_list[h_num])

			# initiate an empty matrix that will be as many rows as the indels file and as many columns as need be
			indels_mat=[] 
			count=0
			
			orig_header_line=''
			indels_reader_list[h_num] = csv.reader(indels_handle_list[h_num], delimiter = '\t')
			for line in indels_reader_list[h_num]:
				count+=1
				if count == 1:
					orig_header_line=line
				
				if count > 1:
					
					if 'A' == 'B': # this is a dummy line to maintain the structure of Gus's original code							
						
						dummy_placeholder = ''
						
					else:
					
					
						indels_mat.append([0]*indels_col_number)	
						for col_index in range(total_cols):
							indels_mat[count-2][col_index]=line[col_index]
						
						c_alt_kmer=''
						for char in indels_mat[count-2][10]:
							c_alt_kmer += comp_dict[char]
	
						c_alt_kmer = c_alt_kmer[::-1]

						kmer_count=0
						for kmernum in range(len(kmerdb_name_list)):
							kmer_count +=1
							if indels_mat[count-2][10] in kmer_topdict['kmer'+str(kmer_count)]:#check to see if the alt kmer is in the dictionary
								indels_mat[count-2][total_cols - 1 + kmer_count]=(kmer_topdict['kmer'+str(kmer_count)][(indels_mat[count-2][10])])
							elif c_alt_kmer in kmer_topdict['kmer'+str(kmer_count)]: #if the alt kmer is not in the dictionary - check to see if the rev comp alt kmer is in the dictionary
								indels_mat[count-2][total_cols - 1 + kmer_count]=(kmer_topdict['kmer'+str(kmer_count)][(c_alt_kmer)])
							else:
								indels_mat[count-2][total_cols - 1 + kmer_count]=0


			if 'A' == 'B': #this is a dummy line to preserve structure of Gus's original codel
				
				dummy_placeholder2 = ''
				
			else:
				
				header_string = ''
				header_string = str(orig_header_line[0])
				for col_header_index in range(1, total_cols):
					header_string += '\t' + str(orig_header_line[col_header_index])
				header_string += kmer_header + '\n'
				outer_low_list[h_num].write(header_string)
				
				for row in indels_mat:
					row_print_string = str(row[0])
					for col_print_index in range(1, indels_col_number):
						row_print_string += '\t' + str(row[col_print_index])
					row_print_string += '\n'
					outer_low_list[h_num].write(row_print_string)
				
	return('\n\ndone\n\n')	
		
if __name__=='__main__':
	print(indels_kmer())