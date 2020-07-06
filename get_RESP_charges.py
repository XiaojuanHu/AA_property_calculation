#!/usr/bin/python
import os, sys, re
import get_RESP_input

#os.system('module load amber/17')

#for folder in os.listdir(os.getcwd()):
#	if not '.' in folder:
#		for little_folder in os.listdir(os.path.join(os.getcwd(), folder)):
		#print folder
		#print os.path.join(path, folder)
		
#get charge
with open(os.path.join(os.getcwd(), 'control.in')) as params:
	file_lines = params.readlines()
	lines = []
	for line in file_lines:
		a = line[0]
		if a != '#':
			lines.append(line)
	
	charge = 'Nan'
	for line in lines:
		if 'charge ' in line:
			rough_charge = line.rsplit(None,1)[-1]
			if '.' in rough_charge:
				charge_1 = rough_charge.split('.')[0]
				charge = int(charge_1)
			else:
				charge = int(rough_charge)
	
	if charge == 'Nan':
		charge = 0.
											
# get pdb
if not os.path.exists(os.path.join(os.getcwd(), 'geometry.pdb')):
	os.system('babel -ifhiaims geometry.in -opdb geometry.pdb')
			
				
			# bunzip potential_esp_1.cube.bz2
		#	if os.path.exists(os.path.join(os.getcwd(), folder, 'potential_esp_1.cube.bz2')):
		#		os.system('cd {} && bunzip2 potential_esp_1.cube.bz2'.format(os.path.join(os.getcwd(), folder)))
		
				
# get esp data file	
if not os.path.exists(os.path.join(os.getcwd(), 'potential_esp_1.cube')):
	sys.exit('== Error: No potential_esp_1.cube found. Exiting now...')
else:
		#		os.system('cd {} && get_RESP_input.py'.format(os.path.join(os.getcwd(), folder)))
	os.system('python ~/scripts/get_RESP_input.py')
		#		os.system('module load amber/17')
	os.system('antechamber -i geometry.pdb -fi pdb -o geometry.ac -fo ac -nc {:d} -pf y -dr no'.format(charge))  # if faulse, try 'antechamber -i geometry.pdb -fi pdb -o geometry.mol2 -fo mol2 -j 5 -at sybyl -dr no'
	os.system('mv ANTECHAMBER_BOND_TYPE.AC0 geometry.ac')
	os.system('respgen -i geometry.ac -o geometry.respin1 -f resp1')
	os.system('resp -O -i geometry.respin1 -o geometry.respout1 -e geometry.esp -t qout_stage1')
	os.system('respgen -i geometry.ac -o geometry.respin2 -f resp2')
	os.system('resp -O -i geometry.respin2 -o geometry.respout2 -e geometry.esp -q qout_stage1 -t qout_stage2')
		
# write charges 
if os.path.exists(os.path.join(os.getcwd(), 'qout_stage2')):
	with open(os.path.join(os.getcwd(), 'qout_stage2')) as structures:
		charges = []
		lines = structures.read()
#			print lines
		for line in lines.split('\n'):
			for i in line.split(None,7):
		#					print i
				charges.append(float(i))
		#			print charges 
						
	with open(os.path.join(os.getcwd(), 'resp.chrg'), 'w') as outputfile:
		for i in charges:
			outputfile.write('{:10.6f}'.format(i) + '\n')
############# clean folder #################################################################################
os.system('bzip2 potential_esp_1.cube')
if not os.path.exists(os.path.join(os.getcwd(), 'resp_files')):
	os.mkdir(os.path.join(os.getcwd(), 'resp_files'))
os.system('mv ATOMTYPE.INF esout geometry.ac geometry.esp geometry.pdb geometry.respin1 geometry.respin2 geometry.respout1 punch qout_stage1 qout_stage2 resp_files')
os.system('tar zcvf resp_files.tar.gz resp_files')
os.system('rm -r resp_files')
	
	
########## collect charges #################################################################################
#if not os.path.exists(os.path.join(os.getcwd(), '..' , 'collect_resp_charges')):
#	os.mkdir(os.path.join(os.getcwd(), '..', 'collect_resp_charges'))
#    
#for folder in os.listdir(os.getcwd()):
#	if not os.path.exists(os.path.join(os.getcwd(), '..', 'collect_resp_charges', folder)):                
#		os.system('mkdir {}'.format(os.path.join(os.getcwd(), '..', 'collect_resp_charges', folder)))
#	os.system('cd {} && cp collect_resp_charges {}'.format(os.path.join(os.getcwd(), folder), os.path.join(os.getcwd(), '..', 'collect_resp_charges', folder)))
				
