#!/usr/bin/python
import os, sys, re

#path = '/home/xhu/Documents/work/DFT/resp/carsten_test/ala/dipeptide/Mg/final'

#unit = 0.52917721
#for folder in os.listdir(os.getcwd()):
#	if not '.' in folder:
#		for little_folder in os.listdir(os.path.join(os.getcwd(), folder)):

# get pdb
if not os.path.exists(os.path.join(os.getcwd(), 'geometry.pdb')):
	os.system('babel -ifhiaims geometry.in -opdb geometry.pdb')
			

# get n_atoms
with open(os.path.join(os.getcwd(), 'geometry.pdb'), 'r') as inputfile:
	file_lines = inputfile.readlines()
	lines = list(file_lines)
			
	n_atoms = 0
	for line in lines:
		coords_found = re.match(r'(\s*(\w+)\s*(\d+)\s*(\w+)\s*(\w+)\s*(\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*(\w+)\s*?)', line)
		if coords_found:
			n_atoms += 1

						
# unzip potential_esp_1.cube.bz2
if os.path.exists(os.path.join(os.getcwd(), 'potential_esp_1.cube.bz2')):
	os.system('bunzip2 potential_esp_1.cube.bz2')
			
sum_esp = 0	
with open(os.path.join(os.getcwd(), 'potential_esp_1.cube'), 'r') as inputfile:
	file_lines = inputfile.readlines()
	lines = list(file_lines)
				
	atom__xyz = []
	lines_position = lines.index(' *****************************\n')
	for line in lines[lines_position+5:lines_position+5+n_atoms]:
		stringa, x, y, z = line.rsplit(None, 3)
		atom__xyz.append([float(x), float(y), float(z)])
	
	# esp : hartree, xyz: bohr	
	esp_xyz = []
	n_esp = 0	
	for line in lines[lines_position+5+n_atoms:]:
		esp_found = re.match(r'(\s*(.?\d+\.\d+\w.?\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+))', line)
		if esp_found:
			potential, x, y, z = line.split(None, 3)
			if not(float(potential) == 0 and float(x) == 0 and float(z) == 0):
				esp_xyz.append([float(potential), float(x), float(y), float(z)])
				n_esp+=1
				sum_esp+=float(potential)
		#	print sum_esp
				
with open(os.path.join(os.getcwd(), 'geometry.esp'), 'w') as outputfile:
	outputfile.write('{:5d}'.format(n_atoms) + '{:5d}'.format(n_esp) +'{:5d}'.format(0) + '\n')
	for i in atom__xyz:
		outputfile.write(' ' + '{:16s}'.format(' ') + '{:16.7E}'.format(i[0]) + '{:16.7E}'.format(i[1]) + '{:16.7E}'.format(i[2]) + '\n')
	for i in esp_xyz:
		outputfile.write(' ' + '{:16.7E}'.format(i[0]*(-1)) + '{:16.7E}'.format(i[1]) + '{:16.7E}'.format(i[2]) + '{:16.7E}'.format(i[3]) + '\n')		
			
		
	
