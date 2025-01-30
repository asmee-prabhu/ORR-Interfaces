#!/usr/bin/env python3

import sys
import numpy as np
import spglib as sym

from ase    	import Atoms
from ase.io 	import read, write
from ase.build 	import sort


# First read input to get supercell 
ifile   = sys.argv[1]
ofile   = sys.argv[2]
slab_sc = read(ifile)

# Check tolerance level required
tol_array = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
for tol in tol_array: 
	spacegroup = sym.get_spacegroup(slab_sc, symprec=tol)
	#print(spacegroup)
	print(f'At a tolerance level of {tol:.2e} A, the space group of the cell is {spacegroup}.')
tol_reqd_f = float(input("Enter tolerance required: "))

# Reduce supercell to primitive cell, creating an ASE Atoms object
unit_cell, scaled_positions, numbers = \
sym.standardize_cell(slab_sc, 
		to_primitive=True, 
		no_idealize=False, 
		symprec=tol_reqd_f,
		)
slab_pc = Atoms(numbers, cell=unit_cell, scaled_positions=scaled_positions)
zpos = slab_pc.get_positions()[:,2]
slab_pc = sort(slab_pc, tags=zpos)
#slab_pc.write(ofile)

# Check if redefinition of lattice is required
if slab_pc.get_chemical_symbols()[-1] == 'Ag' \
or slab_pc.get_chemical_symbols()[-1] == 'Au' \
or slab_pc.get_chemical_symbols()[-1] == 'Pt':
	print('No redefinition of lattice required; primitive cell is written as POSCAR!')
	slab_pc.write(ofile)
else:
	transmat = np.array([[1, 0, 0],
				[0,-1, 0],
				[0, 0,-1]])
	numbers =  numbers
	unit_cell = np.matmul(unit_cell, transmat)
	scaled_positions[:,2] = -1*scaled_positions[:,2]
	slab_pc = Atoms(numbers, cell=unit_cell, scaled_positions=scaled_positions)
	new_unit_cell, new_scaled_positions, new_numbers = \
	sym.standardize_cell(slab_pc, 
			symprec=1e-8,
			)
	slab_pc = Atoms(new_numbers, cell=new_unit_cell, scaled_positions=new_scaled_positions)
	print(f'Lattice redefined; primitive cell is written as {ofile}.')
	slab_pc.write(ofile)
