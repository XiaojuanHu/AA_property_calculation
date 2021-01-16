#!/usr/bin/env python
import os, sys, re, shutil
import matplotlib.pyplot as plt
import fafoom
from fafoom.utilities import atom_masses, VDW_radii
import numpy as np
from itertools import combinations, chain
import math
#from scipy.spatial.distance import pdist, squareform
from collections import OrderedDict
from collections import Counter
import operator, collections

def takeConnectivity(coords, species):
	"""Takes connectivity from atoms and species"""
	
	def testConnect(element, coords):
		def take_vdWs(atoms):
			return 1.1*max([VDW_radii[atoms[0][0]], VDW_radii[atoms[0][1]]])
		def take_length(atoms):
			p0, p1=atoms[1][0], atoms[1][1]
			return math.sqrt((p0[0] - p1[0]) ** 2 + (p0[1] - p1[1]) ** 2 + (p0[2] - p1[2]) ** 2)
		vdWs=[take_vdWs(i) for i in element]
		lengths=[take_length(i) for i in element]
		inxxx=list(combinations(range(len(coords)), 2))
		return [inxxx[i] for i in range(len(vdWs)) if vdWs[i]>=lengths[i]]
	
	atomCobinations =list(combinations(species, 2))
	coordCobinations=list(combinations(coords, 2))
	all_combinations = list(zip(atomCobinations, coordCobinations))
	
	return testConnect(all_combinations, coords)

def constructGraph(conn_list):
	graph = {}
	for i in conn_list:
		if i[0] not in graph:
			graph[i[0]] = [i[1]]
		else:
			graph[i[0]].append(i[1])
	for i in conn_list:
		if i[1] not in graph:
			graph[i[1]] = [i[0]]
		else:
			graph[i[1]].append(i[0])
	return graph
    
def aims2xyz(aims_file):
	coords, species = [], []
	with open(aims_file, 'r') as aims:
		lines = aims.readlines()
		for line in lines:
			atoms = re.match(r'(.*atom\s*?(.?\d+\.\d+)\s+(.?\d+\.\d+)\s+(.?\d+\.\d+)\s+(\w+).*?)', line)
			if atoms:
				coords.append([float(atoms.group(2)), float(atoms.group(3)), float(atoms.group(4))])
				species.append(str(atoms.group(5)))
	return coords, species


#-----------------------------------------------------------------------
# Extract coordinates and atomspecies:
# Output will be array [coord1, coord2, coord3, atommass]
def coords_masses(xyz_file):
	def to_array(ll):
		return([float(ll[1]), float(ll[2]), float(ll[3]), atom_masses[ll[0]]])
	data = [to_array(line.split()) for line in xyz_file]
	return(np.array(data))
	
def coords_masses_light(xyz_file):
	# Takes into account only atoms lighter that Oxygen but not Hydrogens:
	def to_array(ll):
		if 2.0 < atom_masses[ll[0]] < 18.0:
			return([float(ll[1]), float(ll[2]), float(ll[3]), atom_masses[ll[0]]])

	data = [to_array(line.split()) for line in xyz_file if to_array(line.split()) is not None]
	return(np.array(data))	


def coords_masses_light_hs(xyz_file):
	# Takes into account only atoms lighter that Oxygen but not Hydrogens:
	def to_array(ll):
		if atom_masses[ll[0]] < 18.0:
			return([float(ll[1]), float(ll[2]), float(ll[3]), atom_masses[ll[0]]])

	data = [to_array(line.split()) for line in xyz_file if to_array(line.split()) is not None]
	return(np.array(data))	

def coords_masses_heavy(xyz_file):
	# Takes into account only atoms lighter that Oxygen but not Hydrogens:
	def to_array(ll):
		if atom_masses[ll[0]] > 18.0:
			return([float(ll[1]), float(ll[2]), float(ll[3]), atom_masses[ll[0]]])

	data = [to_array(line.split()) for line in xyz_file if to_array(line.split()) is not None]
	return(np.array(data))	

def longest_distance(xyz):
	D = squareform(pdist(xyz))
	N = np.max(D)
	I = np.argmax(D)
	I_row, I_col = np.unravel_index(I, D.shape)
	return I_row, I_col

def dist_between_sets_atoms(xyz1, xyz2):
	dist = [(Pa, Pb, math.sqrt(math.pow(Pa[0]-Pb[0],2) + math.pow(Pa[1]-Pb[1],2) + math.pow(Pa[2]-Pb[2],2)))
		for Pa in xyz1[:, :3] for Pb in xyz2[:, :3]]
	return min(dist, key=lambda x:x[2])[0:2]
	
	
def dist_between_atoms(xyz, atom1, atom2):
	return np.linalg.norm(xyz[atom1][:3] - xyz[atom2][:3])


def COM(coords, species):
	return np.average(coords, axis=0, weights=[atom_masses[i] for i in species])

# ~ def physysorb_chemisorb(coords, species):
    # ~ def separate_mol_surf(coords, species):
        # ~ mol_xyz, surf_xyz=[],[]
        # ~ mol_atoms, surf_atoms=[],[]
        # ~ for i in range(len(species)):
            # ~ if 2.0<atom_masses[species[i]]<18.0:
                # ~ mol_xyz.append(coords[i])
                # ~ mol_atoms.append(species[i])
            # ~ elif atom_masses[species[i]]>18.0:
                # ~ surf_xyz.append(coords[i])
                # ~ surf_atoms.append(species[i])
        # ~ return np.array(mol_xyz), mol_atoms, np.array(surf_xyz), surf_atoms

    
    # ~ conn_list=takeConnectivity(coords, species)
    # ~ graph=constructGraph(conn_list)
    # ~ # Carbons notation:
    # ~ alpha=['C', 'C', 'N', 'H']
    # ~ carboxyl=['O', 'O', 'C']
    # ~ delta=['N', 'C', 'H', 'H']
    # ~ guanidine=['N', 'N', 'N']
    
    # ~ carbons=[alpha, carboxyl, delta, guanidine]
    # ~ target_carbons=[sorted(x) for x in carbons]
    # ~ target_carbons_name=['Alpha', 'Carboxyl', 'Delta', 'Guanidine']
    
    # ~ # Notations
    # ~ Amino_notation={2:'NH2', 3:'NH3'}
    # ~ Carboxyl_notation={0:'COO', 1:'COOH'}
    # ~ backbone=['CH2CH2CH2']
    # ~ Delta_notation={0:'N', 1:'NH'}
    # ~ Guanidine_notation={4:'CNH2NH2', 3:'CNH2NH', 2:'CNHNH'}
    
    # ~ notation_names=[Amino_notation, Carboxyl_notation, Delta_notation, Guanidine_notation]
    # ~ heavy_inds=[i for i in graph if species[i]=='C']
    # ~ summ_hydrogen=7
    # ~ summ_heavy=6
    # ~ for i in graph:
        # ~ con=[species[k] for k in graph[i]]
        # ~ if any([sorted(x)==sorted(con) for x in target_carbons]):
            # ~ ind=target_carbons.index(sorted(con))
            # ~ name=target_carbons_name[ind]  
            # ~ if name=='Alpha':
                # ~ h_ind=graph[i][con.index('H')]
                # ~ if coords[i][-1]>coords[h_ind][-1]:
                    # ~ chiral='down'
                # ~ else:
                    # ~ chiral='up'  
    # ~ mol_xyz, mol_atoms, surf_xyz, surf_atoms = separate_mol_surf(coords, species)
    # ~ shortest_mol_surf=dist_between_sets_atoms(mol_xyz, surf_xyz)
    # ~ ind_mol = mol_xyz.tolist().index(shortest_mol_surf[0].tolist())
    # ~ ind_surf = surf_xyz.tolist().index(shortest_mol_surf[1].tolist())
    # ~ dist = np.linalg.norm(np.array(mol_xyz[ind_mol]) - np.array(surf_xyz[ind_surf]))
    # ~ vdW_shell = (VDW_radii[mol_atoms[ind_mol]] + VDW_radii[surf_atoms[ind_surf]])*0.75
    # ~ dist_above=COM(mol_xyz, mol_atoms)[-1] - max([i[-1] for i in surf_xyz.tolist()])
    # ~ return dist_above, dist, chiral, mol_atoms[ind_mol], surf_atoms[ind_surf]

def extract(lines, name):
	""" Extracts: Notation  Energy COM Shortest Chiral Closest_molatom Closest_surfatom"""   
	COORDS_all, SPECIES_all=[], []
	conn_list, ens, COM, Shortes, chiral, closest_atom, closest_surfatom=[], [], [], [],[], [], []
	#for line in lines:
#		if 'Energy' in line:
#			start=lines.index(line)
#			en=line.split()[1]
#			if  en.endswith('.7'):
#				ens.append(float(en[:-2]))
#			else:
#				ens.append(float(en))
				
	coords, species=[], []
		#atoms=int(lines[start-1])
	for conf_line in lines:
		coords_found = re.match(r'(\s*(\w+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*(.?\d+\.\d+)\s*?)', conf_line)
		if coords_found:
			l=conf_line.split()
			species.append(l[0])
			coords.append(np.array([float(l[1]), float(l[2]), float(l[3])]))
		
	conn_list.append(constructGraph(takeConnectivity(coords, species)))
	COORDS_all.append(coords)
	SPECIES_all.append(species)
	if len(name)>114:
		COM_t, Shortes_t, chiral_t, closest_atom_t, closest_surfatom_t = physysorb_chemisorb(coords, species)
		COM.append(COM_t)
		Shortes.append(Shortes_t)
		chiral.append(chiral_t)
		closest_atom.append(closest_atom_t)
		closest_surfatom.append(closest_surfatom_t)
	if len(name)>411:
		return conn_list, COM, Shortes, chiral, closest_atom, closest_surfatom
	else:
		return conn_list,  COORDS_all, SPECIES_all
			

def extract_proj(lines):
	""" Extracts: Proj1 and Proj2 from lines"""   
	proj1, proj2=[], []
	for line in lines:
		proj1.append(float(line.split()[0]))
		proj2.append(float(line.split()[1]))
	return proj1, proj2

def write_to_cvs(dat, proj_dat, name):
	extracts=['Notation',  'Energy', 'COM', 'Shortest', 'Chiral', 'Closest_molatom', 'Closest_surfatom', 'proj1', 'proj2']
	extracts_mol=['Notation',  'Energy', 'proj1', 'proj2']
	# Writes Information to cvs file
	with open('Hierarchy_{}.cvs'.format(name), 'w') as f:
		if len(name)>4:
			f.write('{}\n'.format(''.join(['{:<40}'.format(z) for z in extracts])))
		else:
			f.write('{}\n'.format(''.join(['{:<40}'.format(z) for z in extracts_mol])))
		# Takes line by line from all elements	
		for k in range(len(dat[0])):
			for i in range(len(dat)):
				f.write('{:<40}'.format(dat[i][k]))
			for p in range(len(proj_dat)):
				f.write('{:<40}'.format(proj_dat[p][k]))
			f.write('\n')
            
def Identify_atomtypes(graph, species):

	""" Return soma atom types with correct names"""
	types_dictionary={}
	residue_dictionary={}
	
	#Lets find Specific Carbons first:
	alpha=sorted(['C', 'C', 'H', 'N'])   
	carbon_ace=sorted(['C', 'H', 'H', 'H'])
	carbon_nme=sorted([ 'H', 'H', 'H', 'N',])
	
	
	target_carbons=[alpha, carbon_ace, carbon_nme]
	target_carbons_name=['Calpha', 'C_ACE', 'C_NME']
	for i in graph:
		# Goes through all connections in graph
		# con is atom species, connected to atom i
		# graph[i] = [connected atoms to atom i in graph]
		# graph is a dictionary {atom:[connections to at1, at2, at3...], atom2:[connections to at1, at2, at3.. ]}
		con=[species[k] for k in graph[i]]
		if sorted(con)==alpha:
			types_dictionary[i]= 'CA'     #'Calpha'
			alpha_ind=i
	
	lista = []			
	for i in graph: 
		con=[species[k] for k in graph[i]]
		if sorted(con)==carbon_ace and alpha_ind not in graph[i]:
			for j in graph[i]:
				if species[j]=='C':
					for k in graph[j]:
						lista.append(species[k])
					if sorted(lista) == ['C', 'N', 'O']:
						types_dictionary[i]='CAY'             # CT_ACE
						#residue_dictionary[i]='ACED'
						C_ACE_ind=i
						#if species[k]=='O': 
						#	types_dictionary[i]='CAY'             # CT_ACE
						#	residue_dictionary[i]='ACED'
						#	C_ACE_ind=i
		elif sorted(con)==carbon_nme and alpha_ind not in graph[i]:
			types_dictionary[i]= 'CAT'       #'C_NME'
			#residue_dictionary[i]='CT3'   #'NME'
			C_NME_ind=i
	
	a = 0	
	for atom in graph[C_ACE_ind]:
		if atom not in types_dictionary.keys():
			if species[atom]=='H':
				a += 1
				types_dictionary[atom]=str('HY' + str(a))        # Hace
				#residue_dictionary[atom]='ACED'
			elif species[atom]=='C':
				types_dictionary[atom]= 'CY'      # 'Cace'
				#residue_dictionary[atom]='ACED'
				ind_Cace=atom
	
	b = 0
	for atom in graph[C_NME_ind]:
		if atom not in types_dictionary.keys():
			if species[atom]=='H':
				b += 1
				types_dictionary[atom]= str('HT' + str(b))          #'Hnme'
				#residue_dictionary[atom]='CT3'
			elif species[atom]=='N':
				types_dictionary[atom]= 'NT'          #'Nnme'
				#residue_dictionary[atom]='CT3'
				ind_Nnme=atom
	
	for atom in graph[ind_Nnme]:
		if atom not in types_dictionary.keys():
			if species[atom]=='H':
				types_dictionary[atom]= 'HNT'          #'HNnme'
				#residue_dictionary[atom]='CT3'
			elif species[atom]=='C':
				types_dictionary[atom]='C'    #'CAA'
	
	for atom in graph[ind_Cace]:
		if atom not in types_dictionary.keys():
			if species[atom]=='O':
				types_dictionary[atom]= 'OY'      #'Oace'
				#residue_dictionary[atom]='ACED'
			elif species[atom]=='N':
				types_dictionary[atom]= 'N'     #'NAA'
	
	ind_NAA= types_dictionary.keys()[types_dictionary.values().index('N')]
	
	for atom in graph[ind_NAA]:
		if atom not in types_dictionary.keys():
			if species[atom]=='H':
				types_dictionary[atom]='HN'   #'HNAA'
			#elif species[atom]=='C':
			#	print ('yes', atom)
			#	types_dictionary[atom]='CA'    #'Calpha'
				
	ind_CAA= types_dictionary.keys()[types_dictionary.values().index('C')]
	for atom in graph[ind_CAA]:
		if atom not in types_dictionary.keys():
			if species[atom]=='O':
				types_dictionary[atom]='O'    #'OAA'
			#elif species[atom]=='C':
			#	types_dictionary[atom]='CA'    #'Calpha'
				
#------------------R group(hisD-HSD)---------------------------------------------------------------------  
	
	ind_Calpha=types_dictionary.keys()[types_dictionary.values().index('CA')]
	for atom in graph[ind_Calpha]:
		if atom not in types_dictionary.keys():
			if species[atom]=='C':
				types_dictionary[atom]='CB'    #'CR'
			elif species[atom]=='H':
				types_dictionary[atom]='HA'    #'Halpha'
	
	ind_CB=types_dictionary.keys()[types_dictionary.values().index('CB')]
	
	if file == 'hisD' or file == 'hisE' or file == 'hisH':		
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]=str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]='CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='C':
					types_dictionary[atom]='CD2'   
				elif species[atom]=='N':
					types_dictionary[atom]='ND1' 
					
		ind_CD2=types_dictionary.keys()[types_dictionary.values().index('CD2')]
		for atom in graph[ind_CD2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]='HD2'
				elif species[atom]=='N':
					types_dictionary[atom]='NE2' 
					
		ind_ND1=types_dictionary.keys()[types_dictionary.values().index('ND1')]
		for atom in graph[ind_ND1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]='HD1'
				elif species[atom]=='C':
					types_dictionary[atom]='CE1' 
					
		ind_CE1=types_dictionary.keys()[types_dictionary.values().index('CE1')]
		for atom in graph[ind_CE1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]='HE1'
				elif species[atom]=='N':
					types_dictionary[atom]='NE2' 
					
		ind_NE2=types_dictionary.keys()[types_dictionary.values().index('NE2')]
		for atom in graph[ind_NE2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]='HE2'
		
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
				
			if count not in residue_dictionary.keys():
				if 'HD1' in types_dictionary.values() and 'HE2' not in types_dictionary.values():
					residue_dictionary[count] = 'HSD'
				elif 'HE2' in types_dictionary.values() and 'HD1' not in types_dictionary.values():
					residue_dictionary[count] = 'HSE'
				elif 'HD1' in types_dictionary.values() and 'HE2' in types_dictionary.values():
					residue_dictionary[count] = 'HSP'
	
	elif file == 'arg' or file == 'argH':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]=str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]='CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		b = 0
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]=str('HG' + str(b))   
				elif species[atom]=='C':
					types_dictionary[atom]='CD'
					
		ind_CD=types_dictionary.keys()[types_dictionary.values().index('CD')]
		c = 0
		for atom in graph[ind_CD]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					c += 1
					types_dictionary[atom]=str('HD' + str(c))   
				elif species[atom]=='N':
					types_dictionary[atom]='NE'
					
		ind_NE=types_dictionary.keys()[types_dictionary.values().index('NE')]
		for atom in graph[ind_NE]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]='HE'  
				elif species[atom]=='C':
					types_dictionary[atom]='CZ'
					
		ind_CZ=types_dictionary.keys()[types_dictionary.values().index('CZ')]
		if file == 'argH':
			d = 0
			for atom in graph[ind_CZ]:
				if atom not in types_dictionary.keys():
					if species[atom]=='N':
						d += 1
						types_dictionary[atom]=str('NH' + str(d)) 
			ind_NH1=types_dictionary.keys()[types_dictionary.values().index('NH1')]
			e = 0
			for atom in graph[ind_NH1]:
				if atom not in types_dictionary.keys():
					if species[atom]=='H':
						e += 1
						types_dictionary[atom]= str('HH1' + str(e))
			ind_NH2=types_dictionary.keys()[types_dictionary.values().index('NH2')]
			f = 0
			for atom in graph[ind_NH2]:
				if atom not in types_dictionary.keys():
					if species[atom]=='H':
						f += 1
						types_dictionary[atom]= str('HH2' + str(f))
		elif file == 'arg':
			for atom in graph[ind_CZ]:
				if atom not in types_dictionary.keys():
					if species[atom]=='N':
						d = 0
						for i in graph[atom]:
							if i not in types_dictionary.keys():
								if species[i]=='H':
									d += 1
						if d == 1:
							types_dictionary[atom]='NH1'
						elif d == 2:
							types_dictionary[atom]='NH2'

					
			ind_NH1=types_dictionary.keys()[types_dictionary.values().index('NH1')]
			for atom in graph[ind_NH1]:
				if atom not in types_dictionary.keys():
					if species[atom]=='H':
						types_dictionary[atom]= 'HH12'
			ind_NH2=types_dictionary.keys()[types_dictionary.values().index('NH2')]
			f = 0
			for atom in graph[ind_NH2]:
				if atom not in types_dictionary.keys():
					if species[atom]=='H':
						f += 1
						types_dictionary[atom]= str('HH2' + str(f))
							
			
				
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'ARG'
					
	elif file == 'asp' or file == 'aspH':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]=str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]='CG'
		
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		if file == 'asp':
			b = 0
			for atom in graph[ind_CG]:
				if atom not in types_dictionary.keys():
					if species[atom]=='O':
						b += 1
						types_dictionary[atom]=str('OD' + str(b))
		elif file == 'aspH':
			for atom in graph[ind_CG]:
				if atom not in types_dictionary.keys():
					if species[atom]=='O':
						for i in graph[atom]:
							if i not in types_dictionary.keys():
								if species[i]=='H':
									types_dictionary[atom] = 'OD2'
			for atom in graph[ind_CG]:
				if atom not in types_dictionary.keys():
					if species[atom]=='O':
						types_dictionary[atom] = 'OD1'
									
			ind_OD2=types_dictionary.keys()[types_dictionary.values().index('OD2')]
			for atom in graph[ind_OD2]:
				if atom not in types_dictionary.keys():
					if species[atom]=='H':
						types_dictionary[atom]= 'HD2'
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'ASP'
	
	elif file == 'ala':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'ALA'
					
	elif file == 'asn':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]='CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='O':
					types_dictionary[atom]= 'OD1'
				elif species[atom]=='N':
					types_dictionary[atom]='ND2'
		
		# cis trans			
		ind_ND2=types_dictionary.keys()[types_dictionary.values().index('ND2')]
		ind_OD1=types_dictionary.keys()[types_dictionary.values().index('OD1')]
		index_distance = {} 
		for atom in graph[ind_ND2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					distance = math.sqrt(   (coords[atom][0] - coords[ind_OD1][0] ) ** 2. \
										+ ( coords[atom][1] - coords[ind_OD1][1] ) ** 2. \
										+ ( coords[atom][2] - coords[ind_OD1][2] ) ** 2. )
					index_distance[atom] = distance
		key_list = []
		for key in index_distance.keys():
			key_list.append(key)
			
		if float(index_distance[key_list[0]]) < float(index_distance[key_list[1]]):
			types_dictionary[key_list[0]] = 'HD21'
			types_dictionary[key_list[1]] = 'HD22'
		else:
			types_dictionary[key_list[1]] = 'HD21'
			types_dictionary[key_list[0]] = 'HD22'
					#b += 1
					#types_dictionary[atom]= str('HD2' + str(b))
			
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'ASN'
					
	elif file == 'cys':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='S':
					types_dictionary[atom]='SG'
		ind_S= types_dictionary.keys()[types_dictionary.values().index('SG')]
		for atom in graph[ind_S]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]='HG1'
		
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'CYS'
					
	elif file == 'gln':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]='CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		b = 0
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HG' + str(b))
				elif species[atom]=='C':
					types_dictionary[atom]='CD'
					
		ind_CD=types_dictionary.keys()[types_dictionary.values().index('CD')]
		for atom in graph[ind_CD]:
			if atom not in types_dictionary.keys():
				if species[atom]=='O':
					types_dictionary[atom]='OE1'
				elif species[atom]=='N':
					types_dictionary[atom]='NE2'
		
		ind_NE2=types_dictionary.keys()[types_dictionary.values().index('NE2')]
		ind_OE1=types_dictionary.keys()[types_dictionary.values().index('OE1')]
		index_distance = {} 
		for atom in graph[ind_NE2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					distance = math.sqrt(   (coords[atom][0] - coords[ind_OE1][0] ) ** 2. \
										+ ( coords[atom][1] - coords[ind_OE1][1] ) ** 2. \
										+ ( coords[atom][2] - coords[ind_OE1][2] ) ** 2. )
					index_distance[atom] = distance
		key_list = []
		for key in index_distance.keys():
			key_list.append(key)
			
		if float(index_distance[key_list[0]]) < float(index_distance[key_list[1]]):
			types_dictionary[key_list[0]] = 'HE21'
			types_dictionary[key_list[1]] = 'HE22'
		else:
			types_dictionary[key_list[1]] = 'HE21'
			types_dictionary[key_list[0]] = 'HE22'
		
					
		#ind_NE2=types_dictionary.keys()[types_dictionary.values().index('NE2')]
		#b = 0
		#for atom in graph[ind_NE2]:
		#	if atom not in types_dictionary.keys():
		#		if species[atom]=='H':
		#			b += 1
		#			types_dictionary[atom]= str('HE2' + str(b))
		
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'GLN'
					
	elif file == 'glu' or file == 'gluH':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]='CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		b = 0
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HG' + str(b))
				elif species[atom]=='C':
					types_dictionary[atom]='CD'
					
		ind_CD=types_dictionary.keys()[types_dictionary.values().index('CD')]
		if file == 'glu':
			b = 0
			for atom in graph[ind_CD]:
				if atom not in types_dictionary.keys():
					if species[atom]=='O':
						b += 1
						types_dictionary[atom]= str('OE' + str(b))
		elif file == 'gluH':
			for atom in graph[ind_CD]:
				if atom not in types_dictionary.keys():
					if species[atom]=='O':
						for i in graph[atom]:
							if i not in types_dictionary.keys():
								if species[i]=='H':
									types_dictionary[atom] = 'OE2'
			for atom in graph[ind_CD]:
				if atom not in types_dictionary.keys():
					if species[atom]=='O':
						types_dictionary[atom] = 'OE1'
										
			ind_OE2=types_dictionary.keys()[types_dictionary.values().index('OE2')]
			for atom in graph[ind_OE2]:
				if atom not in types_dictionary.keys():
					if species[atom]=='H':
						types_dictionary[atom]= 'HE2'
		
		
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'GLU'
					
	elif file == 'ile':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HB'
				elif species[atom]=='C':
					for i in graph[atom]:
						if i not in types_dictionary.keys():
							if species[i]=='C':
								types_dictionary[atom]= 'CG1'
								
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='C':
					types_dictionary[atom]= 'CG2'
					
		ind_CG2=types_dictionary.keys()[types_dictionary.values().index('CG2')]
		b = 0
		for atom in graph[ind_CG2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HG2' + str(b))
		
		ind_CG1=types_dictionary.keys()[types_dictionary.values().index('CG1')]
		c = 0
		for atom in graph[ind_CG1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					c += 1
					types_dictionary[atom]= str('HG1' + str(c))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CD'
					
		ind_CD=types_dictionary.keys()[types_dictionary.values().index('CD')]
		c = 0
		for atom in graph[ind_CD]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					c += 1
					types_dictionary[atom]= str('HD' + str(c))
			
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'ILE'
	
	if file == 'leu':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		b = 0
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HG'
				elif species[atom]=='C':
					b += 1
					types_dictionary[atom]= str('CD' + str(b))
					
		ind_CD1=types_dictionary.keys()[types_dictionary.values().index('CD1')]
		b = 0
		for atom in graph[ind_CD1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HD1' + str(b))
		
		ind_CD2=types_dictionary.keys()[types_dictionary.values().index('CD2')]
		b = 0
		for atom in graph[ind_CD2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HD2' + str(b))
		
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'LEU'
	
	elif file == 'lys' or file == 'lysH':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		b = 0
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HG' + str(b))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CD'
					
		ind_CD=types_dictionary.keys()[types_dictionary.values().index('CD')]
		b = 0
		for atom in graph[ind_CD]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HD' + str(b))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CE'
					
		ind_CE=types_dictionary.keys()[types_dictionary.values().index('CE')]
		b = 0
		for atom in graph[ind_CE]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HE' + str(b))
				elif species[atom]=='N':
					types_dictionary[atom]= 'NZ'
					
		ind_NZ=types_dictionary.keys()[types_dictionary.values().index('NZ')]
		b = 0
		for atom in graph[ind_NZ]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HZ' + str(b))
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'LYS'
					
	elif file == 'met':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		b = 0
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HG' + str(b))
				elif species[atom]=='S':
					types_dictionary[atom]= 'SD'
					
		ind_SD=types_dictionary.keys()[types_dictionary.values().index('SD')]
		for atom in graph[ind_SD]:
			if atom not in types_dictionary.keys():
				if species[atom]=='C':
					types_dictionary[atom]= 'CE'
					
		ind_CE=types_dictionary.keys()[types_dictionary.values().index('CE')]
		b = 0
		for atom in graph[ind_CE]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					b += 1
					types_dictionary[atom]= str('HE' + str(b))
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'MET'
					
	elif file == 'phe':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CG'
				
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		b = 0
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='C':
					b += 1
					types_dictionary[atom]= str('CD' + str(b))
		
		ind_CD1=types_dictionary.keys()[types_dictionary.values().index('CD1')]
		for atom in graph[ind_CD1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HD1'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CE1'
		
		ind_CD2=types_dictionary.keys()[types_dictionary.values().index('CD2')]
		for atom in graph[ind_CD2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HD2'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CE2'
					
		ind_CE1=types_dictionary.keys()[types_dictionary.values().index('CE1')]
		for atom in graph[ind_CE1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HE1'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CZ'
		
		ind_CE2=types_dictionary.keys()[types_dictionary.values().index('CE2')]
		for atom in graph[ind_CE2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HE2'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CZ'
					
		ind_CZ=types_dictionary.keys()[types_dictionary.values().index('CZ')]
		for atom in graph[ind_CZ]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HZ'
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'PHE'
	
	elif file == 'pro':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CG'
				
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		a = 0
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HG' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CD'
					
		ind_CD=types_dictionary.keys()[types_dictionary.values().index('CD')]
		a = 0
		for atom in graph[ind_CD]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HD' + str(a))
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'PRO'
	
	elif file == 'ser':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='O':
					types_dictionary[atom]= 'OG'
				
		ind_OG=types_dictionary.keys()[types_dictionary.values().index('OG')]
		for atom in graph[ind_OG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HG1'
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'SER'
					
	elif file == 'thr':
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HB'
				elif species[atom]=='O':
					types_dictionary[atom]= 'OG1'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CG2'
					
		ind_OG1=types_dictionary.keys()[types_dictionary.values().index('OG1')]
		for atom in graph[ind_OG1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HG1'
					
		ind_CG2=types_dictionary.keys()[types_dictionary.values().index('CG2')]
		a = 0
		for atom in graph[ind_CG2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HG2' + str(a))
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'THR'
	
	elif file == 'trp':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='C':
					for i in graph[atom]:
						if i not in types_dictionary.keys():
							if species[i]=='H':
								types_dictionary[atom]= 'CD1'
								
		ind_CD1=types_dictionary.keys()[types_dictionary.values().index('CD1')]
		for atom in graph[ind_CD1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HD1'
				elif species[atom]=='N':
					types_dictionary[atom]= 'NE1'
					
		ind_NE1=types_dictionary.keys()[types_dictionary.values().index('NE1')]
		for atom in graph[ind_NE1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HE1'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CE2'
					
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='C':
					types_dictionary[atom]= 'CD2'
					
		ind_CD2=types_dictionary.keys()[types_dictionary.values().index('CD2')]
		for atom in graph[ind_CD2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='C':
					types_dictionary[atom]= 'CE3'
					
		ind_CE3=types_dictionary.keys()[types_dictionary.values().index('CE3')]
		for atom in graph[ind_CE3]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HE3'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CZ3'
					
		ind_CZ3=types_dictionary.keys()[types_dictionary.values().index('CZ3')]
		for atom in graph[ind_CZ3]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HZ3'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CH2'
					
		ind_CH2=types_dictionary.keys()[types_dictionary.values().index('CH2')]
		for atom in graph[ind_CH2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HH2'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CZ2'
					
		ind_CZ2=types_dictionary.keys()[types_dictionary.values().index('CZ2')]
		for atom in graph[ind_CZ2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HZ2'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CE2'
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'TRP'
					
	elif file == 'tyr':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HB' + str(a))
				elif species[atom]=='C':
					types_dictionary[atom]= 'CG'
					
		ind_CG=types_dictionary.keys()[types_dictionary.values().index('CG')]
		a = 0
		for atom in graph[ind_CG]:
			if atom not in types_dictionary.keys():
				if species[atom]=='C':
					a += 1
					types_dictionary[atom]= str('CD' + str(a))
					
		ind_CD1=types_dictionary.keys()[types_dictionary.values().index('CD1')]
		for atom in graph[ind_CD1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HD1'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CE1'
					
		ind_CD2=types_dictionary.keys()[types_dictionary.values().index('CD2')]
		for atom in graph[ind_CD2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HD2'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CE2'
					
		ind_CE1=types_dictionary.keys()[types_dictionary.values().index('CE1')]
		for atom in graph[ind_CE1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HE1'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CZ'
					
		ind_CE2=types_dictionary.keys()[types_dictionary.values().index('CE2')]
		for atom in graph[ind_CE2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HE2'
				elif species[atom]=='C':
					types_dictionary[atom]= 'CZ'
					
		ind_CZ=types_dictionary.keys()[types_dictionary.values().index('CZ')]
		for atom in graph[ind_CZ]:
			if atom not in types_dictionary.keys():
				if species[atom]=='O':
					types_dictionary[atom]= 'OH'
					
		ind_OH=types_dictionary.keys()[types_dictionary.values().index('OH')]
		for atom in graph[ind_OH]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					types_dictionary[atom]= 'HH'
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'TYR'
					
	elif file == 'val':
		a = 0
		for atom in graph[ind_CB]:
			if atom not in types_dictionary.keys():
				if species[atom]=='C':
					a += 1
					types_dictionary[atom]= str('CG' + str(a))
				elif species[atom]=='H':
					types_dictionary[atom]= 'HB'
					
		ind_CG1=types_dictionary.keys()[types_dictionary.values().index('CG1')]
		a = 0
		for atom in graph[ind_CG1]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HG1' + str(a))
					
		ind_CG2=types_dictionary.keys()[types_dictionary.values().index('CG2')]
		a = 0
		for atom in graph[ind_CG2]:
			if atom not in types_dictionary.keys():
				if species[atom]=='H':
					a += 1
					types_dictionary[atom]= str('HG2' + str(a))
					
		for count, item in enumerate(species):
			if item == 'Ba':
				types_dictionary[count] = 'BA'
				residue_dictionary[count] = 'BA'
			
			if count not in residue_dictionary.keys():
					residue_dictionary[count] = 'VAL'
	
	
	
	
	
	if len(types_dictionary.keys())==len(species):
		print 'Done!'
		#print  (types_dictionary)
		#print residue_dictionary
	return types_dictionary, residue_dictionary


if __name__=='__main__':
	"""Make description what script does"""
	file = str(sys.argv[1])
	filename = file + '.xyz'
	#path='/home/damaksimovda/Temp/Connectivity_convert/'
	#coords, species= aims2xyz(os.path.join(path, 'geometry.in'))
	#trajectory = open(os.path.join(path, filename)).readlines()
	#path='/home/xhu/scripts/conn_convert/'
	coords, species= aims2xyz(os.path.join(os.getcwd(),'geometry.in'))
	trajectory = open(os.path.join(os.getcwd(),'coords.xyz')).readlines()
	
	def species_identify(conn, species):
		species_related={}
		for i in conn:
			species_related[species[i]]=[species[k] for k in conn[i]] 
		return species_related
			
	conn_dict,  COORDS, SPECIES=extract(trajectory, file)
	for i in range(len(conn_dict)):
		types_dictionary, residue_dictionary = Identify_atomtypes(conn_dict[i], SPECIES[i])
	#residues = []	
	#for value in residue_dictionary.values():
	#	if value not in residues:
	#		residues.append(value)
			
	#atom_group = {}
	#for key in residue_dictionary.keys():
	#	if residue_dictionary[key] == residues[0]:
	#		atom_group[key] = 1
	#	elif residue_dictionary[key] == residues[1]:
	#		atom_group[key] = 2
	#	elif residue_dictionary[key] == residues[2]:
	#		atom_group[key] = 3
	#	elif residue_dictionary[key] == residues[3]:
	#		atom_group[key] = 4
		
	#print types_dictionary
	#print residue_dictionary
	#print atom_group
	
	#generate pdb file for CHARMM
	with open((os.path.join(os.getcwd(), 'residue.pdb')), 'w') as pdbfile:
		pdbfile.write('CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1' + '\n')
		
		for key in types_dictionary.keys()[:-1]:  # If there is no cation in the system, then for key in types_dictionary.keys():
			pdbfile.write('{:6s}{:>5d}'.format('ATOM', int(key) + 1)+' ' +'{:4s}'.format(types_dictionary[key])+' ' + '{:<4s}'.format(residue_dictionary[key])+ 'P' + '{:4d}'.format(int(1))\
			 + '    '+'{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}'.format(float(coords[key][0]), float(coords[key][1]), float(coords[key][2]), float(0.00),float(0.00))+'      '\
			 +'PROA' + '\n' )
		pdbfile.write('END'+'\n'+'\n')
	
	# If there is cation in the system
	with open((os.path.join(os.getcwd(), 'ba.pdb')), 'w') as pdbfile:
		pdbfile.write('CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1' + '\n')
		
		pdbfile.write('{:6s}{:>5d}'.format('ATOM', int(types_dictionary.keys()[-1]) + 1)+' ' +'{:4s}'.format('BA')+' ' + '{:<4s}'.format('BA')+ 'P' + '{:4d}'.format(int(2))\
			+ '    '+'{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}'.format(float(coords[types_dictionary.keys()[-1]][0]), float(coords[types_dictionary.keys()[-1]][1]), float(coords[types_dictionary.keys()[-1]][2]), float(0.00),float(0.00))+'      '\
			+'HETA' + '\n' )
		pdbfile.write('END'+'\n'+'\n')	
		
	#generate pdb file for openmm
	#with open((os.path.join(os.getcwd(), 'openmm.pdb')), 'w') as pdbfile:
	#	pdbfile.write('REMARK    GENERATED BY xhu' + '\n')
	#	pdbfile.write('HEADER    Dipeptide in gas' + '\n')
	#	pdbfile.write('MODEL        0' + '\n')
	#	
	#	for key in types_dictionary.keys():
	#		pdbfile.write('{:6s}{:>5d}'.format('ATOM', int(key) + 1)+ ' ' +'{:4s}'.format(types_dictionary[key])+' ' + '{:3s}'.format(residue_dictionary[key]) + '  '+ \
	#		'{:4d}'.format(residue_index[key]) + '    '+'{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}'.format(float(coords[key][0]), float(coords[key][1]), float(coords[key][2]), float(1.00),float(0.00))\
	#		+'          '+'{:2s}'.format(species[key])+'\n')
	#
	#	pdbfile.write('END'+'\n'+'\n')
		
	sys.exit(0)
	conn_dict=constructGraph(takeConnectivity(coords, species))
	# ~ print conn_dict
	# ~ for i in conn_dict:
		# ~ for k in conn_dict[i]:
			# ~ print i, k
	# ~ print protonation_state(coords, species)
	
	
	
	
	# Path to folder with data
	data_path='/home/damaksimovda/Nextcloud/Papers/Arg_on_surfaces/structure_space_all/cutoff_3.0/Sketchmaps_projections/mds/data/'
	
	# Names of systems:
	# ~ names=['Arg', 'ArgH', 'ArgCu', 'ArgAg', 'ArgAu',  'ArgHCu', 'ArgHAg', 'ArgHAu']
	names=['ArgH']
	

	# Create dictionary {name:lines of sorted xyz_file}
	dictionary=OrderedDict()
	for i in names:
		dictionary[i]=open(os.path.join(data_path, '{}_sorted.xyz'.format(i))).readlines()
	
	# Creates list with len()=len(names) 
	# each element contains all the information 
	# about all structures in xyz trajectory file
	data=[list(extract(dictionary[i], i)) for i in dictionary]

    #Extract projections
	proj_data=[list(extract_proj(open(os.path.join(data_path, '{}.proj'.format(i))).readlines())) for i in dictionary]
	# Writes all extracted data to cvs files
	for i in range(len(data)):
		write_to_cvs(data[i], proj_data[i], names[i])
