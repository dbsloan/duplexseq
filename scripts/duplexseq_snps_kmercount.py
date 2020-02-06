#!/usr/bin/env python3

import sys
import csv
import re
import os
import copy

def snps_kmer():

	'''
	This function uses a kmer database(s) and a set of SNPs from our duplexseq processing pipeline to report kmer abundances
	
	It is called with:
	NameOfScript
		followed by these 5 arguments (all separated with spaces)
		1. input file with path
		2. kmer database path 
		3. list of kmer databases separated by commas (no spaces). For example: symb_shotgun_k39max10000_short.txt,symb_shotgun_k39max10000_short2.txt
		   if user only has one kmer database, that is fine, put spaces before and after and don't use any commas
		4. output file and path
		5. Total number of columns in input file
	
	
#example call: duplexseq_snps_kmercount.py  SNPs_file.txt  /home/KMC_databases/  my_kmer_db.txt  SNPs_file.kmer.txt  13
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
		kmer_header = kmer_header + '\tRefcount_' + kmer_db_file + '\tAltcount_' + kmer_db_file
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
		
			
	snps_col_number =total_cols+(2*(len(kmerdb_name_list)))
	
		
	###
	###
	###
	'''
	loop through the snps files, find kmers, and write to ouput
	'''
	###
	###
	###
	
	
	
	#create lists of file handles and output files for each user specified rep
	snps_counter=0
	
	#there is no longer a need for a loop here, but I am keeping the structure of gus's original code intact
	for dummy_var in [1]:
		
		snps_file_list=[]
		snps_handle_list=[]
		snps_reader_list=[]
		outer_low_list=[]
	
		snps_counter+=1
		snps_file_list.append(inputfile)
		snps_handle_list.append('file' + str(snps_counter))
		snps_reader_list.append(('reader' + str(snps_counter)))
		outer_low_list.append(outputfile)
		
		
		for h_num in range(len(snps_handle_list)):
			try:
				snps_handle_list[h_num]=open(snps_file_list[h_num])		
				outer_low_list[h_num]=open(outer_low_list[h_num], 'w')

			except:
				return ('unable to open/create the file: ' + snps_file_list[h_num] + ' or ' + outer_low_list[h_num])

			# initiate an empty matrix that will be as many rows as the snps file and as many columns as need be
			snps_mat=[] 
			count=0
			
			orig_header_line=''
			snps_reader_list[h_num] = csv.reader(snps_handle_list[h_num], delimiter = '\t')
			for line in snps_reader_list[h_num]:
				count+=1
				if count == 1:
					orig_header_line=line
				
				if count > 1:
					
					if 'A' == 'B': # this is a dummy line to maintain the structure of Gus's original code							
						
						dummy_placeholder = ''
						
					else: # it is a snps file
					
					
						snps_mat.append([0]*snps_col_number)	
						for col_index in range(total_cols):
							snps_mat[count-2][col_index]=line[col_index]
						
						c_ref_kmer=''
						for char in snps_mat[count-2][10]:
							c_ref_kmer += comp_dict[char]
						
						c_ref_kmer = c_ref_kmer[::-1]
						
						c_alt_kmer=''
						for char in snps_mat[count-2][11]:
							c_alt_kmer += comp_dict[char]
	
						c_alt_kmer = c_alt_kmer[::-1]

						kmer_count=0
						for kmernum in range(len(kmerdb_name_list)):
							kmer_count +=1
							if snps_mat[count-2][10] in kmer_topdict['kmer'+str(kmer_count)]:#check to see if the ref kmer is in the dictionary
								snps_mat[count-2][total_cols - 2 + (2*(kmer_count))]=(kmer_topdict['kmer'+str(kmer_count)][(snps_mat[count-2][10])])
							elif c_ref_kmer in kmer_topdict['kmer'+str(kmer_count)]: #if the ref kmer is not in the dictionary - check to see if the rev comp ref kmer is in the dictionary
								snps_mat[count-2][total_cols - 2 + (2*(kmer_count))]=(kmer_topdict['kmer'+str(kmer_count)][(c_ref_kmer)])
							else:
								snps_mat[count-2][total_cols - 2 + (2*(kmer_count))]=0

							if snps_mat[count-2][11] in kmer_topdict['kmer'+str(kmer_count)]:#check to see if the alt kmer is in the dictionary
								snps_mat[count-2][total_cols - 1 + (2*(kmer_count))]=(kmer_topdict['kmer'+str(kmer_count)][(snps_mat[count-2][11])])
							elif c_alt_kmer in kmer_topdict['kmer'+str(kmer_count)]: #if the alt kmer is not in the dictionary - check to see if the rev comp alt kmer is in the dictionary
								snps_mat[count-2][total_cols - 1 + (2*(kmer_count))]=(kmer_topdict['kmer'+str(kmer_count)][(c_alt_kmer)])
							else:
								snps_mat[count-2][total_cols - 1 + (2*(kmer_count))]=0


			if 'A' == 'B': #this is a dummy line to preserve structure of Gus's original codel
				
				dummy_placeholder2 = ''
				
			else:
				
				header_string = ''
				header_string = str(orig_header_line[0])
				for col_header_index in range(1, total_cols):
					header_string += '\t' + str(orig_header_line[col_header_index])
				header_string += kmer_header + '\n'
				outer_low_list[h_num].write(header_string)
				
				for row in snps_mat:
					row_print_string = str(row[0])
					for col_print_index in range(1, snps_col_number):
						row_print_string += '\t' + str(row[col_print_index])
					row_print_string += '\n'
					outer_low_list[h_num].write(row_print_string)
				
	return('\n\ndone\n\n')	
		
if __name__=='__main__':
	print(snps_kmer())