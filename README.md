# AA_property_calculation
Scripts for property calculation of amino acid data set. 

## get_cube_params.py
Get parameters for the generation of total density cube file and electrostatic potential cube file in FHI-aims, and write the relevent comments into the parameter file of FHI-aims: 'control.in'.

## bader.py
Derive bader charges from the results of bader charge analysis code (http://theory.cm.utexas.edu/henkelman/code/bader/).

## get_RESP_input.py
Get antechamber input files: geometry.pdb(generate from geometry.in by openbabel), geometry.esp(electrostatic potential file in antechamber input format).

## get_RESP.charges.py
Script to derive resp charges. Antechamber and get_RESP_input.py are called.

## data_collection.sh
Collect energy, geometry, partial charges(bader, hirshfeld, resp) information. Side files: coord.xyz(geometry in xyz format), hirsh.chrg(hirshfeld charges).
