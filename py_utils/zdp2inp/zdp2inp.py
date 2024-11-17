############################################################################################
#                                                                                          #
# Author: Praise Noah Johnson                                                              #
#                                                                                          #
# Description: This code generates a .dat file for 0D plasma-assisted combustion data,     #
#              extracts the raw data for PAC from various data files, and also generates   #
#              a chem.inp file with all the reactions (combustion + dummy plasma reactions)#
#                                                                                          #
############################################################################################
import cantera as ct
import numpy as np
import sys
import csv
import re
import os
import math as mt
import os.path
from os import path
from shutil import copyfile





# Extracts the plasma reaction from kinet.inp and creates a chem.inp file with all the reactions.
# The rates of plasma reactions are zeros.
def extract_plasma_reactions(plasma_reactions, inp_name):
	
	# Opening a new .inp file to write
	if path.exists(inp_name):
	
		print('Opening mech_name to write plasma reactions...')
		inp_file = open('chem.inp','w')
		
		section = ''
		
		# Opening the existing combustion mechanism to copy the stuff into new file
		for line in open(inp_name, 'r'):
			if re.match('^REACTIONS', line):
				section = 'reactions'
			
			# Copying the combustion reactions
			if section != 'reactions':
				inp_file.write(line)
			elif section == 'reactions':
				if re.match('^END', line):
					continue
				else:
					inp_file.write(line)
	
	else:
	
		print('Specified mech_name file not found...')
		sys.exit()
	
	
	if path.exists(plasma_reactions):
		print('Start reading plasma_reactions to extract plasma reactions...')
	else:
		print ('Specified plasma_reactions file not found...')
		sys.exit()
	
	
	plasma_reaction_list_cti = []
	plasma_reaction_list_inp = []
	section = ''
	previous_line = ''
	count = 0
		
	
	# Opening the kinet.inp file to read the plasma reactions
	for line in open(plasma_reactions, 'r'):
		
		# Skipping lines starting with #, $, or blank lines
		if re.match('^#', line):
			continue
		elif line == '^\n':
			continue
		elif line == '':
			continue
		elif re.match('^\$', line):
			continue
		elif re.match('^\s+$', line):
			continue
		
		# Clearing newlines, inline comments, and endline spaces.
		line = re.sub('\n','',line)
		line = re.sub('!.+','',line)
		line = re.sub('\s*$','',line)
		
		
		# Ckeck if reactions has started
		if re.match('^REACTIONS', line):
			section = 'reactions'
			continue
		
		
		# Converting the plasma reactions into chemkin format dummy reaction lines
		if section == 'reactions':
			
			split_line = []
			
			# See if current line starts like @B = ... If so, print previous reaction line for corresponding times.
			if re.match('^\s*@B\s*=\s*', line):
				
				split_line = re.split('\s+', re.sub('\s*@B\s*=\s*','',line))
				
				for item in split_line:
				
					count += 1
					print (str(count)+'\t'+re.sub('@B',item,previous_line))
					plasma_reaction_list_cti.append('reaction(\'' +  previous_line.replace('@B', item).replace('\'', '\\\'') + '\', [0.0, 0.0, 0.0])')
					plasma_reaction_list_inp.append(previous_line.replace('@B', item) + '\t\t0.0 0.0 0.0')
			
			# See if current line starts like @M = ... If so, print previous reaction line for corresponding times.
			elif re.match('^\s*@M\s*=\s*', line):
				
				split_line = re.split('\s+', re.sub('\s*@M\s*=\s*','',line))
				
				for item in split_line:
				
					count += 1
					print (str(count)+'\t'+re.sub('@M',item,previous_line))
					plasma_reaction_list_cti.append('reaction(\'' +  previous_line.replace('@M', item).replace('\'', '\\\'') + '\', [0.0, 0.0, 0.0])')
					plasma_reaction_list_inp.append(previous_line.replace('@M', item) + '\t\t0.0 0.0 0.0')
			
			# See if current line starts like @R = ...
			elif re.match('^\s*@R\s*=\s*', line):
				
				# See if previous line starts like @B = ... or @M = ... If not, print previous reaction line for corresponding times.
				if not (re.match('^\s*@B\s*=\s*', previous_line) or re.match('^\s*@M\s*=\s*', previous_line)):
					split_line = re.split('\s+', re.sub('\s*@R\s*=\s*','',line))
					
					for item in split_line:
					
						count += 1
						print (str(count)+'\t'+previous_line)
						print ('\t\tDup')
						plasma_reaction_list_cti.append('reaction(\'' +  previous_line.replace('\'', '\\\'') + '\', [0.0, 0.0, 0.0])')
						plasma_reaction_list_inp.append(previous_line + '\t\t0.0 0.0 0.0\n\tDup')
			
			# Provision to not miss the last reaction
			elif re.match('^END', line):
				
				count += 1
				print (str(count)+'\t'+previous_line)
				plasma_reaction_list_cti.append('reaction(\'' +  previous_line.replace('\'', '\\\'') + '\', [0.0, 0.0, 0.0])')
				plasma_reaction_list_inp.append(previous_line + '\t\t0.0 0.0 0.0')
			
			# Provision to not miss the first reaction
			else:
				
				if not (re.match('^\s*@B\s*=\s*', previous_line) or re.match('^\s*@M\s*=\s*', previous_line) or re.match('^\s*@R\s*=\s*', previous_line)):
					
					if previous_line != '':
					
						count += 1
						print (str(count)+'\t'+previous_line)
						plasma_reaction_list_cti.append('reaction(\'' +  previous_line.replace('\'', '\\\'') + '\', [0.0, 0.0, 0.0])')
						plasma_reaction_list_inp.append(previous_line + '\t\t0.0 0.0 0.0')
			
			
			previous_line = re.sub('\n','',line)
	
	
	# Ckeck if number of plasma reactions matches with the plasma reaction rate data
	#if len(net_plasma_rr) != len(plasma_reaction_list_cti):
	#	print len(net_plasma_rr), "\t", len(plasma_reaction_list_cti)
	#
	#	print 'Mismatch in number of plasma reactions and plasma reaction rate data. Check if plasma reaction file is formatted correctly...'
	#	sys.exit()
		
	print ('\nWriting the following lines to the chem.inp mechanism file...')
	
	
	# Generating a chem.inp file with all the reactions
	for i in range(len(plasma_reaction_list_cti)):
	
		#print plasma_reaction_list_cti[i]
		print (str(i+1) + '\t' + plasma_reaction_list_inp[i])
		#mech_file.write('\n\n# Plasma Reaction ' + str(i+1) + '\n')
		#mech_file.write(plasma_reaction_list_cti[i])
		#mech_file.write('\n# Dummy reaction written for GPSA')
		inp_file.write('\n\n! Plasma Reaction ' + str(i+1) + '\n')
		inp_file.write(plasma_reaction_list_inp[i])
		inp_file.write('\n! Dummy reaction written for GPSA')
	
	inp_file.write('\n\nEND')
	
	
	#mech_file.close()
	inp_file.close()
	print ('\t\t\t\t\t\t\t-----------------------> Done')
	






extract_plasma_reactions("kinet.inp", "chem1.inp")
