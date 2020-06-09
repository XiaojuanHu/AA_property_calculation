#!/usr/bin/python
import os, sys, re
import pandas as pd
import numpy
import math

# read geometry.in file
with open(os.path.join(os.getcwd(), 'geometry.in')) as structures:
	lines = structures.readlines()

atom_x = []
atom_y = []
atom_z = []
for line in lines:
	coords_found = re.match(r'(\s*(\w+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*(\w+)\s*?)', line)
	if coords_found:
		stringa, x, y, z, ele = line.split(None)
		atom_x.append(float(x))
		atom_y.append(float(y))
		atom_z.append(float(z))
		
# 14 bohr -> angstrom
len_x = (max(atom_x)- min(atom_x)) + 14*0.529177249
len_y = (max(atom_y)- min(atom_y)) + 14*0.529177249
len_z = (max(atom_z)- min(atom_z)) + 14*0.529177249

# get steps along each directions
x_steps = math.floor(len_x/0.05)
y_steps = math.floor(len_y/0.05)
z_steps = math.floor(len_z/0.05)
#print(x_steps, y_steps, z_steps)

# write cube parameter sentences to control.in
new_line_x = ' cube edge ' + str(x_steps) + '  0.0500    0.0000    0.0000'
new_line_y = ' cube edge ' + str(y_steps) + '  0.0000    0.0500    0.0000'
new_line_z = ' cube edge ' + str(z_steps) + '  0.0000    0.0000    0.0500'

# total density parameters
os.system('sed -i "/^#  max_relaxation_steps/a\ cube_content_unit bohr" control.in')
os.system('sed -i "/ cube_content_unit bohr/i\#############################################" control.in')
os.system('sed -i "/ cube_content_unit bohr/i\ output hirshfeld" control.in')
os.system('sed -i "/ cube_content_unit bohr/i\#" control.in')
os.system('sed -i "/ cube_content_unit bohr/i\ output cube total_density" control.in')
os.system('sed -i "/ cube_content_unit bohr/i\{:s}" control.in'.format(new_line_x))
os.system('sed -i "/ cube_content_unit bohr/i\{:s}" control.in'.format(new_line_y))
os.system('sed -i "/ cube_content_unit bohr/i\{:s}" control.in'.format(new_line_z))

# esp parameters
os.system('sed -i "/ cube_content_unit bohr/a\ esp output_cube 6" control.in')
os.system('sed -i "/ esp output_cube 6/i\#" control.in')
os.system('sed -i "/ esp output_cube 6/i\ output esp" control.in')
os.system('sed -i "/ esp output_cube 6/i\ esp radius 1.4 2.0" control.in')
os.system('sed -i "/ esp output_cube 6/i\ esp equal_grid_n 35 35 35" control.in')
