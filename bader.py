#!/usr/bin/python
import os, sys, re
import pandas as pd
import numpy

# This script is to derive bader charges from the results of bader charge analysis code.

# get electron number
ACF_file = open('ACF.dat', 'r')
ACF = list(ACF_file)
ACF_file.close()

electron_num = []
for line in ACF:
	charge_found = re.match(r'(\s*(\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*?)', line)
	if charge_found:
		n, x, y, z, chrg, min_dist, vol = line.split(None)
		electron_num.append(float(chrg))
		
# write out electron number
with open('bdr.chrg', 'w') as bdr:
	for i in electron_num:
		bdr.write('{:10.6f}\n'.format(i))
		

# get nuclear charge
in_file = open('output', 'r')
file_lines = in_file.readlines()
lines = list(file_lines)
in_file.close()

for line in lines:
        if '  The structure contains' in line:
                n_atoms = int(line.split(None, 4)[-2])
        
lines_position = lines.index('  Input geometry:\n')
        
for line in lines[lines_position+4:lines_position+4+n_atoms]:
	stringa, n_atom, stringb, symbol, x, y, z = line.split(None, 6)
	elements.append(symbol)



nuc = []
for i in elements:
	if i == 'C':
		nu = 6
	elif i == 'O':
		nu = 8
	elif i == 'H':
		nu = 1
	elif i == 'N':
		nu = 7
	elif i == 'S':
		nu = 16
	elif i == 'Na':
		nu = 11
	elif i == 'Mg':
		nu = 12
	elif i == 'P':
		nu = 15
	elif i == 'Ca':
		nu = 20
	elif i == 'Zn':
		nu = 30
		
	nuc.append(nu)

with open('nuclear.chrg', 'w') as out:
	for i in nuc:
		out.write('{:10.6f}\n'.format(i))
		
# get bader charge
bader = []
for i in range(len(nuc)):
	ba = nuc[i] - float(electron_num[i])
	bader.append(ba)
	
with open('bader.chrg', 'w') as out:
	for i in bader:
		out.write('{:10.6f}\n'.format(i))
		
		
		
		
		
		
		
		
		
		
		
		
		
		
