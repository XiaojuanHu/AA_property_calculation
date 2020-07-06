#!/usr/bin/bash

# get energies.dat
echo "# Calculated with PBE+evdW (tier 2 basis)" > energies.dat
echo "# conf of   Energy/ev" >> energies.dat

for i in conf_*; do 
	E1=`grep 'Total energy of the DFT' $i/output | head -n 1 | awk '{print $12}'`
	echo $i $E1 >> energies.dat
	done

# get geometry.ext
for i in conf_*; do 
	cd $i
	
	grep atom geometry.in | awk '{print $5,$2,$3,$4}' > coords.xyz
	
	## Hirshfeld charges
	grep ' Hirshfeld charge        :' output | awk '{print $5}' > hirsh.chrg
	
	echo "#atom: x : y : z : bader : hirsf. : resp" >> geometry.ext
	paste coords.xyz bader.chrg hirsh.chrg resp.chrg >> geometry.ext
	
	done
