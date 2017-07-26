#!/usr/bin/python
### Python Library for VASP
### INPUT: CONTCAR, OUTCAR, DOSCAR 
### Last update: 12 February 2014
from copy import deepcopy
import os
import numpy as np
path = os.getcwd().split('/')[-1]

###CONTCAR
def get_total_number_of_atoms():
	with open('CONTCAR', 'r') as cont:
		for i in range(6): tmp = cont.readline()
		num_atoms = 0
		a = cont.readline().split()
		for n in range(len(a)):
			num_atoms += int(a[n])
		return num_atoms

def get_number_of_each_atom():
	with open('CONTCAR', 'r') as cont:
		for i in range(6): tmp = cont.readline()
		return cont.readline().split()

def get_sort_of_atoms():
	with open('CONTCAR', 'r') as cont:
		for i in range(5): tmp = cont.readline()
		return cont.readline().split()

def get_atom_coordinates(filename):

    if isinstance(filename, str):
        f = open(filename)
    else:  # Assume it's a file-like object
        f = filename

    line1 = f.readline()
    lattice_constant = float(f.readline().split()[0])
    # Now the lattice vectors
    a = []
    for ii in range(3):
        s = f.readline().split()
        floatvect = float(s[0]), float(s[1]), float(s[2])
        a.append(floatvect)

    basis_vectors = np.array(a) * lattice_constant

    atom_symbols = []
    numofatoms = f.readline().split()

    vasp5 = False
    try:
        int(numofatoms[0])
    except ValueError:
        vasp5 = True
        atomtypes = numofatoms
        numofatoms = f.readline().split()
        numbers = [int(a) for a in numofatoms]

    # check for comments in numofatoms line and get rid of them if necessary
    commentcheck = np.array(['!' in s for s in numofatoms])
    if commentcheck.any():
        # only keep the elements up to the first including a '!':
        numofatoms = numofatoms[:np.arange(len(numofatoms))[commentcheck][0]]

    if not vasp5:
        atomtypes = line1.split()

        numsyms = len(numofatoms)
        if len(atomtypes) < numsyms:
            # First line in POSCAR/CONTCAR didn't contain enough symbols.

            # Sometimes the first line in POSCAR/CONTCAR is of the form
            # "CoP3_In-3.pos". Check for this case and extract atom types
            if len(atomtypes) == 1 and '_' in atomtypes[0]:
                atomtypes = get_atomtypes_from_formula(atomtypes[0])
            else:
                atomtypes = atomtypes_outpot(f.name, numsyms)
        else:
            try:
                for atype in atomtypes[:numsyms]:
                    if atype not in chemical_symbols:
                        raise KeyError
            except KeyError:
                atomtypes = atomtypes_outpot(f.name, numsyms)

    for i, num in enumerate(numofatoms):
        numofatoms[i] = int(num)
        [atom_symbols.append(atomtypes[i]) for na in range(numofatoms[i])]

    # Check if Selective dynamics is switched on
    sdyn = f.readline()
    selective_dynamics = sdyn[0].lower() == 's'

    # Check if atom coordinates are cartesian or direct
    if selective_dynamics:
        ac_type = f.readline()
    else:
        ac_type = sdyn
    cartesian = ac_type[0].lower() == 'c' or ac_type[0].lower() == 'k'
    direct = ac_type[0].lower() == 'd'
    tot_natoms = sum(numofatoms)
    atoms_pos = np.empty((tot_natoms, 3))
    if selective_dynamics:
        selective_flags = np.empty((tot_natoms, 3), dtype=bool)
    for atom in range(tot_natoms):
        ac = f.readline().split()
        atoms_pos[atom] = (float(ac[0]), float(ac[1]), float(ac[2]))
        if selective_dynamics:
            curflag = []
            for flag in ac[3:6]:
                curflag.append(flag == 'F')
            selective_flags[atom] = curflag

    # Done with all reading
    if isinstance(filename, str):
        f.close()
    if direct:
        atoms_pos = [ np.dot(a, basis_vectors) for a in atoms_pos]

    return(atoms_pos, basis_vectors, numbers)

###OUTCAR
def get_keyword_outcar(keyword):
	with open('OUTCAR', 'r') as out:
		for line in out:
			if keyword in line:
				return line.split()[2]

###DOSCAR
def split_DOSCAR(number_of_each_atom, sort_of_atoms, ISPIN, FORBITAL):
	##Skipping unnecessary part in DOSCAR
	doscar = open('DOSCAR')
	for i in range(5): tmp=doscar.readline()
	num_line = int(doscar.readline().split()[2])
	doscar.close

	doscar = open("DOSCAR", "r")
	for i in range(6): tmp = doscar.readline()

	#TDOS, considering ISPIN
	tdos = open("tDOS", "w")
	for i in range(num_line): tdos.write(doscar.readline().strip() + "\n")
	doscar.readline()

	#PDOS, considering either ISPIN and FORBITAL
	for n in range(len(number_of_each_atom)):
		atom = sort_of_atoms[n]
		print "found atom is " + atom
		for i in range(int(number_of_each_atom[n])):
			file_name = "pDOS_" + atom + "_" + str(i + 1)
			pdos = open(file_name, "w")
			for k in range(num_line): pdos.write(doscar.readline().strip() + "\n")
			doscar.readline()
	################################################
	############# writing itx for Igor #############
	################################################
	#tDOS.itx
	tdos_src = open("tDOS", "r")
	tdos_itx = open("tDOS.itx", "w")
	tdos_lines = tdos_src.readlines()
	tdos_itx.write("IGOR\n")
	if ISPIN == True:
		tdos_itx.write("WAVES/D " + "E_tdos_" + path + " TDOS_up_" + path + " TDOS_dw_" + path + " accum_up_" + path + " accum_dw_" + path)	
	else:
		tdos_itx.write("WAVES/D " + "E_tdos_" + path + " TDOS_" + path + " accum_" + path)		
	tdos_itx.write("\nBEGIN\n")
	for line in tdos_lines: tdos_itx.write(line)
	tdos_itx.write("END\n")
	if ISPIN == True:
		tdos_itx.write("X Display " + "TDOS_up_" + path + " vs " + "E_tdos_" + path + " as " + '"' + path + '"' + "\n")
	else:
		tdos_itx.write("X Display " + "TDOS_" + path + " vs " + "E_tdos_" + path + " as " + '"' + path + '"' + "\n")
	tdos_itx.write("""X ModifyGraph width=340.157,height=340.157
X ModifyGraph marker=19
X ModifyGraph lSize=1.5
X ModifyGraph tick=2
X ModifyGraph mirror=1
X ModifyGraph zero(bottom)=8
X ModifyGraph fSize=20
X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10
X ModifyGraph standoff=0
X ModifyGraph axThick=1.5
X ModifyGraph axisOnTop=1
X Label left "\\Z20 Density-of-states (arb. unit)"
X Label bottom "\\Z20 Energy (eV)"
"""
)

def layer_DOS(numbers, n_bins, a_DOS):
    '''
        numbers: the number of each of the atom types, list of integers
        n_bins: Number of bins that you wish to use for layers (int)
	a_DOS: The list from the atomic_DOS function

    Returns:
	s_l_DOS: A list of layer and species resolved DOS. The format is [species, layer, [DOS]]
    '''

    d_length = len(a_DOS[0][2])
    o_length = len(a_DOS[0][2][0])
    n_spec = len(numbers)

    s_l_DOS = []
    for s in range(len(numbers)):
        l_DOS = []
        for n in range(n_bins):
	    tmp = np.zeros(shape=(len(a_DOS[0][2]),len(a_DOS[0][2][0])))
            for dos in a_DOS:
	        tdos = np.asarray(dos[2])
  	        if dos[1] == n and dos[0] == s:
		    for i in range(len(tmp)):
		        tmp[i,0] = float(tdos[i,0])
		        for j in range(1,len(tmp[0])):
			    tmp[i,j] = tmp[i,j] + float(tdos[i,j])
            l_DOS.append(tmp)
        s_l_DOS.append(l_DOS)
 

    return(s_l_DOS) 
	    


def atomic_DOS(number_of_each_atom, n_bins, direction, positions, lattice):
    '''
	number_of_each_atom: the number of each of the atom types, list of integers
	n_bins: Number of bins that you wish to use for layers (int)
 	direction: The direction along which you are binning (string)
   	positions: The atomic coordinates
    '''

    if direction == 'x': 
        coord = 0
        length = np.sqrt(sum([x**2 for x in lattice[:,0]]))
    if direction == 'y': 
        coord = 1
        length = np.sqrt(sum([x**2 for x in lattice[:,1]]))
    if direction == 'z': 
        coord = 2
        length = np.sqrt(sum([x**2 for x in lattice[:,2]]))

    bin_width = length / n_bins
    

    doscar = open('DOSCAR')
    for i in range(5): tmp=doscar.readline()
    num_line = int(doscar.readline().split()[2])
    doscar.close

    doscar = open("DOSCAR", "r")
    for i in range(6): tmp = doscar.readline()
    for i in range(num_line): tmp = doscar.readline()
    doscar.readline()

    atoms_dos = []
    for n in range(len(number_of_each_atom)):
        for i in range(int(number_of_each_atom[n])):
	    tmp_dos = []
            for j in range(num_line): tmp_dos.append(doscar.readline().split())
	    tmp = doscar.readline()
	    for k in range(n_bins):
		if positions[i][coord] >= bin_width * k and positions[i][coord] < bin_width * ( k + 1 ):
	            atom_dos = [n, k, tmp_dos]
            atoms_dos.append(atom_dos)

    return(atoms_dos) 
    

def sum_DOSCAR(atom, range_min, range_max):
	##make list of pdos file names
	file_list = []
	for a in range(range_min, range_max + 1):
		file_list.append("pDOS_" + atom + "_" + str(a))

	##make empty list for new list, named sum_dos
	pdos_first = [i.strip().split() for i in open("pDOS_" + atom + "_" + str(range_min)).readlines()]
	sum_dos = deepcopy(pdos_first)
	
	for i in range(len(sum_dos)):
		for k in range(len(sum_dos[0])):
			sum_dos[i][k] = 0

	##process of making summation of lists		
	print "Sum of ..."
	for file_name in file_list:
		print file_name + " ", 
		pdos = [i.strip().split() for i in open(file_name).readlines()]
	##convert all numbers from string to float
		for a in range(len(pdos)):
			for b in range(len(pdos[a])):
				pdos[a][b] = float(pdos[a][b])
	##sum of each element in all lists
		for n, line in enumerate(pdos):
			sum_dos[n] = map(sum, zip(line, sum_dos[n]))
	##correction energy values
	for i in range(len(pdos_first)):
		sum_dos[i][0] = pdos_first[i][0]

	##Output
	with open("sum_dos_" + atom + str(range_min) +"to" + str(range_max), "w") as output:
		output.writelines(' '.join(str(j) for j in i ) + "\n" for i in sum_dos)
	
	################################################
	############# writing itx for Igor #############
	################################################
	
