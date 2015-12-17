import os, copy, sys, random, cPickle, math
sys.path.append("/fs/home/jms875/Library")
import utils, filetypes, lammps
import numpy as np
import time

coul_on = True

def center_of_geometry(atoms):
	return utils.Struct(  x = sum(a.x for a in atoms)/len(atoms),
		y = sum(a.y for a in atoms)/len(atoms),
		z = sum(a.z for a in atoms)/len(atoms)   )

def radius_of_gyration(atoms):
		center = center_of_geometry(atoms)
		return math.sqrt( 1.0/len(atoms) * sum( [utils.dist_squared(a,center) for a in atoms] ) )

def build_dot(N, spoc, monomer, atoms, bonds, angles, dihedrals, random_seed=1):
	random.seed(random_seed)
	pb_type = monomer.atoms[0].type
	s_type = monomer.atoms[1].type
	Q = 2.968
	L = int( math.ceil( (2*N)**0.333) )
	center = utils.Struct(x=0., y=0., z=0.)
	core_atoms = []
	if L*Q<30: #octahedron
		L+=6
		for zi in range(L):
			#L2 = int(( math.sqrt(3)*( L/2-abs(zi-L/2) ) ))
			L2 = (L/2-abs(zi-L/2))*2
			for xi in range( L2 ):
				for yi in range( L2 ):	
					x = (xi-L2*0.5+0.5)*Q
					y = (yi-L2*0.5+0.5)*Q
					z = (zi-L*0.5+0.5)*Q
					core_atoms.append( utils.Struct(x=x, y=y, z=z, element='Pb' if (xi+yi+zi)%2 else 'S', type=pb_type if (xi+yi+zi)%2 else s_type) )
	else: #cube
		L+=2
		for xi in range(L):
			for yi in range(L):
				for zi in range(L):
					x = (xi-L*0.5+0.5)*Q
					y = (yi-L*0.5+0.5)*Q
					z = (zi-L*0.5+0.5)*Q
					core_atoms.append( utils.Struct(x=x, y=y, z=z, element='Pb' if (xi+yi+zi)%2 else 'S', type=pb_type if (xi+yi+zi)%2 else s_type) )
	
	#filetypes.write_xyz('out', core_atoms)
	#sys.exit()
	
	for a in core_atoms:
		a.dist = dist=utils.dist(a,center)
		a.neighbors = []
		for b in core_atoms:
			if a is not b and utils.dist_squared(a,b)<Q**2+0.1:
				a.neighbors.append(b)
	
	core_atoms.sort(key=lambda a:a.dist)
	s_atoms = [a for a in core_atoms if a.element=='S'][:N]
	pb_atoms = {}
	for s in s_atoms:
		for a in s.neighbors:
			if a.type is pb_type:
				pb_atoms[a] = True
	
	pb_atoms = pb_atoms.keys()
	pb_atoms.sort(key=lambda a:a.dist)
	#excess_pb_atoms = pb_atoms[N:] #tends to pick all from one side
	
	closest_dist_to_remove = pb_atoms[N].dist
	
	core_pb_atoms = [a for a in pb_atoms if a.dist<closest_dist_to_remove-0.01]
	marginal_pb_atoms = [a for a in pb_atoms if ( abs(a.dist-closest_dist_to_remove)<0.01 and a not in core_pb_atoms) ]
	random.shuffle(marginal_pb_atoms)
	core_pb_atoms = core_pb_atoms + marginal_pb_atoms[ : N-len(core_pb_atoms) ]
	excess_pb_atoms = [a for a in pb_atoms if a not in core_pb_atoms]

	for a in s_atoms + core_pb_atoms:
		atoms.append(a)
	
	import numpy

	def add_spoc_at_pb(pb):
		spoc.add_to(0.0,0.0,0.0, atoms, bonds, angles, dihedrals)
	
		direction = [pb.x, pb.y, pb.z]
		direction /= numpy.linalg.norm(direction)
		
		dot = numpy.dot(direction,[1., 0., 0.])
		theta = math.acos(dot)
		cross = numpy.cross(direction,[1., 0., 0.])
		cross /= -numpy.linalg.norm(cross)
		w = numpy.array(cross)
		
		for a in atoms[-len(spoc.atoms):]:
			v = numpy.array([a.x, a.y, a.z])
			a.x, a.y, a.z = v*math.cos(theta) + numpy.cross(w,v)*math.sin(theta) + w*numpy.dot(w,v)*(1-math.cos(theta)) #Rodrigues' rotation formula
			a.x += pb.x
			a.y += pb.y
			a.z += pb.z
	
	for pb in excess_pb_atoms:
		add_spoc_at_pb(pb)
	
	#add more Pb complexes as needed
	
	new_pb = 0
	while new_pb < N+3:
		x,y,z = [ Q*( 3+ (2*N)**0.333)/2 ]*3
		M = utils.rand_rotation()
		x,y,z = utils.matvec(M, [x,y,z])
		pb = utils.Struct(x=x, y=y, z=z, element='Pb', type=pb_type)
		too_close = False
		for a in atoms:
			if utils.dist_squared(a,pb) < (Q+2)**2: too_close = True; break
		if not too_close: add_spoc_at_pb(pb); new_pb+=1
	
	#filetypes.write_xyz('out', atoms)
	#sys.exit()

c18_ligands = False

PbS = 909
PbO = 907
S = 908
OO = 214
OH = 376
HO = 377
COO = 213
if c18_ligands:
	spoc_name = 'pb-c18-oh'
	CCO = 216
else:
	spoc_name = 'pb-e-oh'
	CCO = 215
HC = 85

def write_input_header(f, run_name, atom_types, read_restart=False):
	type_indices = dict( [(t.index,i+1) for i,t in enumerate(atom_types)] )
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] )
	pb_types = [i+1 for i,t in enumerate(atom_types) if t.element==82]
	s_types = [i+1 for i,t in enumerate(atom_types) if t.element==16]
	o_types = [i+1 for i,t in enumerate(atom_types) if t.element==8]
	c_types = [i+1 for i,t in enumerate(atom_types) if t.element==6]
	h_types = [i+1 for i,t in enumerate(atom_types) if t.element==1]
	if not read_restart:
		f.write('''units		real
atom_style	full
variable pairstyle string '''+('lj/cut/coul/long' if coul_on else 'lj/cut')+'''
pair_style hybrid/overlay lj/cut'''+('/coul/long 10.0' if coul_on else ' 10.0')+''' morse 10.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary p p p
read_data	'''+run_name+'''.data''')

	f.write('''
	kspace_style pppm 1.0e-3
'''+( '\n'.join([ ('#%d:\t%d\t%s' % (i+1,t.index,t.element)) for i,t in enumerate(atom_types) ]) )+'''

pair_coeff * * morse 0.0 1.0 1.0
pair_coeff * * ${pairstyle}	0.0	1.0\n''')

	#set normal LJ from OPLS
	f.write(('\n'.join(["pair_coeff %d %d ${pairstyle}	%f %f" % (atom_type_numbers[t], atom_type_numbers[t], t.vdw_e, t.vdw_r) for i,t in enumerate(atom_types) if i+1 not in pb_types+s_types])) + '\n')

	#for stability, prevent crashes:
	for pair in [ (OO,CCO), (OH,CCO) ]:
		f.write('pair_coeff	%d	%d	${pairstyle}	0.01	4.5\n' % tuple(sorted([type_indices[i] for i in pair ])) )
	for pair in [ (OH,HC), (OO,HC)]:
		f.write('pair_coeff	%d	%d	${pairstyle}	0.01	3.5\n' % tuple(sorted([type_indices[i] for i in pair ])) )
	#for pair in [ (S,S), (S,OO), (S,OH), (OH,OH) ]:
	#	f.write('pair_coeff	%d	%d	${pairstyle}	0.01	3.5\n' % tuple(sorted([type_indices[i] for i in pair ])) )
	for i in c_types: #repel Pb and S from C
		for j in pb_types+s_types:
			f.write('pair_coeff	%d	%d	${pairstyle}	0.01	4.0\n' % tuple(sorted([i, j])) )
	for i in h_types: #repel Pb and S from H
		for j in pb_types+s_types:
			f.write('pair_coeff	%d	%d	${pairstyle}	0.01	3.5\n' % tuple(sorted([i, j])) )
	
	for types in params: #LJ params - overwrite any earlier ones
		if type(types)==type(1) or len(params[types])!=2: continue
		if types[0] in type_indices and types[1] in type_indices:
			f.write('pair_coeff %d %d ${pairstyle} %f %f \n' % tuple( sorted([type_indices[i] for i in types]) + params[types] ) )
	
	for types in params: #Morse params
		if type(types)==type(1) or len(params[types])!=3: continue
		if types[0] in type_indices and types[1] in type_indices:
			f.write('pair_coeff %d %d morse %f %f %f\n' % tuple( sorted([type_indices[i] for i in types]) + params[types] ) )
			f.write('pair_coeff %d %d ${pairstyle} 0.0 1.0 \n' % tuple( sorted([type_indices[i] for i in types]) ) ) #set corresponding LJ to zero

def run_job(run_name, on_queue=True):
	if on_queue:
		f = open(run_name+'.nbs', 'w')
		f.write('''#!/bin/bash
##NBS-nproc: 1
##NBS-queue: "batch"

/fs/home/jms875/build/lammps/lammps-9Dec2014/src/lmp_serial -in %s.in -log %s.log &> /dev/null
''' % (run_name, run_name) )
		f.close()
		os.system("jsub "+run_name+".nbs")
	else: #run locally
		os.system('/fs/home/jms875/build/lammps/lammps-9Dec2014/src/lmp_serial -in %s.in -log %s.log' % (run_name,run_name) )

def save_to_file( atoms, bonds, angles, dihedrals, atom_types, run_name ):
	#prevent loops in stored object by replacing references with interger indices
	for x in bonds+angles+dihedrals:
		x.atoms = [a.index for a in x.atoms]
	for a in atoms:
		try:
			a.bonded = [b.index for b in a.bonded]
		except: pass
		a.type = a.type.index
		a.neighbors = None
		for x in a.__dict__: #turn things with bad types into floats
			t = type(a.__dict__[x])
			if t not in [type(1), type(1.0), type(''), type(None)]:
				try:
					a.__dict__[x] = float( a.__dict__[x] )
				except:
					pass
	cPickle.dump( (atoms, bonds, angles, dihedrals, atom_types), open(run_name+'.pickle', 'wb') )
	#reset
	for x in bonds+angles+dihedrals:
		x.atoms = [atoms[ii-1] for ii in x.atoms]
	for a in atoms:
		try:
			a.bonded = [atoms[ii-1] for ii in a.bonded]
		except: pass
		a.type = [t for t in atom_types if t.index==a.type][0]

def load_from_file(run_name):
	atoms, bonds, angles, dihedrals, atom_types = cPickle.load( open('lammps/'+run_name+'.pickle', 'rb') )

	#reload references from indices
	for x in bonds+angles+dihedrals:
		x.atoms = [atoms[ii-1] for ii in x.atoms]
	for a in atoms:
		try:
			a.bonded = [atoms[ii-1] for ii in a.bonded]
		except: pass
		a.type = [t for t in atom_types if t.index==a.type][0]
	
	#reload cartesian coordinates from xyz file
	xyz_lines = open('lammps/'+run_name+'.xyz').readlines()[-len(atoms):]
	xyz = [ line.split()[1:] for line in xyz_lines ]
	for i in range(len(xyz)):
		atoms[i].x, atoms[i].y, atoms[i].z = [float(s) for s in xyz[i]]
	#center core atoms
	center = center_of_geometry( [a for a in atoms if a.element=='S'] )
	for a in atoms:
		a.x -= center.x; a.y -= center.y; a.z -= center.z
	return atoms, bonds, angles, dihedrals, atom_types


def remove_unattached_ligands(atoms, bonds, angles, dihedrals, spoc):
	sulfurs = [a for a in atoms if a.element=='S']
	remove_count = 0
	atoms_to_remove = {}
	for i,a in enumerate(atoms):
		if a.element=='Pb':
			nearest_sulfur = min([utils.dist_squared(a,b) for b in sulfurs])
			if nearest_sulfur > 5**2:
				#remove whole spoc, based on size of spoc and where pb is within it (always first?)
				remove_count = len(spoc.atoms)
		if remove_count > 0:
			remove_count -= 1
			atoms_to_remove[a] = True
	atoms = [a for a in atoms if a not in atoms_to_remove]
	bonds = [b for b in bonds if b.atoms[0] not in atoms_to_remove and b.atoms[1] not in atoms_to_remove ]
	angles = [a for a in angles if a.atoms[0] not in atoms_to_remove and a.atoms[1] not in atoms_to_remove and a.atoms[2] not in atoms_to_remove ]
	dihedrals = [d for d in dihedrals if d.atoms[0] not in atoms_to_remove and d.atoms[1] not in atoms_to_remove and d.atoms[2] not in atoms_to_remove and d.atoms[3] not in atoms_to_remove ]
	return atoms, bonds, angles, dihedrals
	

def react_dots(params, dot_size1, dot_size2, prior_run_name1, prior_run_name2, random_seed=1, run_name_prefix='react_', on_queue=True, method='metadynamics'):
	random_s = str(random_seed)
	random.seed(random_seed)
	run_name = run_name_prefix+str(dot_size1)+'_'+str(dot_size2)
	
	spoc = utils.Molecule('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate_morse/'+spoc_name+'.arc')
	
	atoms1, bonds1, angles1, dihedrals1 = cPickle.load( open('lammps/'+prior_run_name1+'.pickle', 'rb') )
	atoms1, bonds1, angles1, dihedrals1 = remove_unattached_ligands(atoms1, bonds1, angles1, dihedrals1, spoc)
	atoms2, bonds2, angles2, dihedrals2 = cPickle.load( open('lammps/'+prior_run_name2+'.pickle', 'rb') )
	atoms2, bonds2, angles2, dihedrals2 = remove_unattached_ligands(atoms2, bonds2, angles2, dihedrals2, spoc)
	
	atoms1_xyz = filetypes.parse_xyz(prior_run_name1+'.xyz')
	atoms2_xyz = filetypes.parse_xyz(prior_run_name2+'.xyz')
	goal_atoms = filetypes.parse_xyz('dot' + str(dot_size1+dot_size2)+'.xyz')
	
	# set starting coordinates
	for a,b in zip(atoms1,atoms1_xyz):
		a.x = b.x
		a.y = b.y
		a.z = b.z
	
	for a,b in zip(atoms2,atoms2_xyz):
		a.x = b.x
		a.y = b.y
		a.z = b.z
		
	atoms = atoms1+atoms2
	bonds = bonds1+bonds2
	angles = angles1+angles2
	dihedrals = dihedrals1+dihedrals2
	# make sure atom types are the same in atom1 and atom2
	atom_types_by_type_index = dict( [(t.type.index,t.type) for t in atoms] )
	for i,a in enumerate(atoms):
		a.type = atom_types_by_type_index[ a.type.index ]
		a.index = i+1
		if coul_on and a.type.index in params:
			a.charge = params[a.type.index]
		else:
			a.charge = a.type.charge
	
	atom_types = dict( [(t.type,True) for t in atoms] ).keys()
	atom_types.sort(key=lambda t:-t.element-t.index*0.00001)
	
	offset = radius_of_gyration( [a for a in atoms1 if (a.type.notes == 'PbS Nanocrystal Pb' or a.element == 'S')] ) + radius_of_gyration( [a for a in atoms2 if (a.type.notes == 'PbS Nanocrystal Pb' or a.element == 'S')] ) + 10.0
	offset_vector = utils.matvec( utils.rand_rotation(), [offset, offset, offset] )
	rot = utils.rand_rotation()
	for a in atoms2: #offset by random orientation vector
		a.x, a.y, a.z = utils.matvec( rot, [a.x, a.y, a.z] )
		a.x += offset_vector[0]
		a.y += offset_vector[1]
		a.z += offset_vector[2]
	box_size = [ 15.0+max([a.x for a in atoms])-min([a.x for a in atoms]),
				 15.0+max([a.y for a in atoms])-min([a.y for a in atoms]),
				 15.0+max([a.z for a in atoms])-min([a.z for a in atoms]) ]
	
	#goal_atoms = cPickle.load( open('lammps/dot'+str(dot_size1+dot_size2)+'.pickle', 'rb') )[0]
	
	#filetypes.write_xyz('out_meta', atoms)
	#exit()
	
	os.chdir('lammps')
	save_to_file( atoms, bonds, angles, dihedrals, atom_types, run_name )
	lammps.write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name, atom_types=atom_types, pair_coeffs_included=False)
	os.system('cp ../'+sys.argv[0]+' '+run_name+'.py')
	
	starting_rxn_coord = radius_of_gyration( [a for a in atoms if (a.type.notes == 'PbS Nanocrystal Pb' or a.element == 'S')] )
	print len([a for a in atoms if (a.type.notes == 'PbS Nanocrystal Pb' or a.element == 'S')])
	goal_rxn_coord = radius_of_gyration( goal_atoms[0:2*(dot_size1+dot_size2)] )
	
	colvars_file = open(run_name+'.colvars', 'w')
	colvars_file.write('''
colvarsTrajFrequency 10000
colvarsRestartFrequency 100000

colvar {
  name gyr
  width 0.1 #size of bins, in Angstroms
  upperBoundary '''+str( starting_rxn_coord )+''' #two separate dots
  lowerBoundary '''+str( goal_rxn_coord )+''' #one unified dot
  lowerwallconstant 10.0
  upperwallconstant 10.0
  gyration {
	atoms {
	  atomNumbers '''+(' '.join([str(a.index) for a in atoms if (a.type.notes == 'PbS Nanocrystal Pb' or a.element == 'S')]))+'''
	}
  }
}
''')
	if method.lower()=='abf':
		colvars_file.write('''
abf {
  colvars gyr
  hideJacobian # when using distance-based colvar: makes pmf flat at infinity
  fullSamples 1000
}
''')
	elif method.lower()=='metadynamics':
		colvars_file.write('''
metadynamics {
  colvars gyr
  hillWeight 0.1
  newHillFrequency 1000
}
''')
	else: raise Exception('Invalid free energy method "%s"' % method)
	colvars_file.close()
	
	f = open(run_name+'.in', 'w')
	write_input_header(f, run_name, atom_types, read_restart=False)
	f.write('''
dump	1 all xyz 10000 '''+run_name+'''.xyz
minimize 0.0 1.0e-8 10000 1000
neigh_modify check yes every 1 delay 0
thermo_style custom pe temp
thermo 10000
restart 100000 '''+run_name+'''.restart1 '''+run_name+'''.restart2
fix motion all nve
fix implicit_solvent all langevin 400.0 400.0 100.0 '''+random_s+''' zero yes gjf yes
fix		col all colvars '''+run_name+'''.colvars output '''+run_name+''' seed '''+random_s+''' tstat implicit_solvent
velocity	all create 400.0 '''+random_s+''' mom yes
timestep	2.0
thermo_style	custom step temp etotal pe ke epair ebond f_col tpcpu
run '''+str(int(1e7))+'''
write_restart '''+run_name+'''.restart
''')
	f.close()
	#run_job(run_name, on_queue)
	os.chdir('..')

params = {
(PbS,PbS):[3.147584,0.4223836,7.5],
(S,S):[0.8750986,0.6291512,6.6],
(PbS,S):[33.09405,2.086741,3.204742],
(PbS,OH):[11.20058,1.850411,2.331682],
(PbS,PbO):[0.01,1.01162,7.091745],
(S,OH):[0.001,1.10159,6.542388],
(S,PbO):[9.636918,1.739724,2.891541],
(PbS,OO):[1.465489,1.251277,2.716893],
(S,OO):[0.01289595,0.7470482,7.5],
(PbO,OH):[25.0864,2.308311,2.130475],
(PbO,PbO):[0.007665716,0.93812,7.5],
(OH,OH):[0.001,2.4,4.045742],
(PbO,OO):[8.694149,2.316947,2.134054],
(OO,OH):[0.001,2.298382,4.3651],
(OO,OO):[0.5,2.4,2],
(OO,COO):[0.8195998,1.810373,3.635776],
(OO,CCO):[0.001505031,1.438961,5.8],
(OH,COO):[0.1,2.5,1.947915],
(OH,CCO):[0.001,2.4,2.2],
(PbS,COO):[0.5,1.975883,3.549393],
(PbO,COO):[0.01,2.388686,4.389351],
(S,COO):[3.85033,2.150051,3.909868],
(PbS,CCO):[3.101,2.2,3.328943],
(PbO,CCO):[3.635663,2.1,3.372447],
(S,CCO):[0.22539,2.409581,3.6],
(COO,COO):[0.5,2.151459,3.260709],
(PbS,HO):[0.0765805,1.22558,4.5],
(PbO,HO):[0.007476037,1.570133,4.9],
(S,HO):[0.0100603,2.2,3.675594],
(OH,HO):[2.368954,1.625292,1.8],
(OO,HO):[0.3797819,1.701958,2.2],
(COO,HO):[0.08145238,2.401465,3.589952],
(COO,CCO):[3.804006,2.5,2.74505],

PbS	: 1.124218,
S	: -1.124218,
PbO	:0.8656056,
OH	:-0.7269209,
OO	:-0.6981355,
HO	:0.370847,
CCO	: -0.28,
HC	: 0.06				
}

if c18_ligands:
	params[COO] = -(params[PbO] + params[OH] + 2*params[OO] + params[HO] + params[CCO] + 2*params[HC])
else:
	params[COO] = -(params[PbO] + params[OH] + 2*params[OO] + params[HO] + params[CCO] + 3*params[HC]) #one more H to account for in ethanoate ligands

utils.Molecule.set_params('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate_morse/oplsaa.prm')

def react_all_dots():
	for step in range(1):
		for sum_size in range(2,max_dot_size+1):
			for dot_size1 in range(1,sum_size/2+1):
				dot_size2 = sum_size-dot_size1
				react_dots(params, dot_size1, dot_size2, 'c18_'+str(dot_size1), 'c18_'+str(dot_size2), random_seed=step+1, run_name_prefix='react_%d__' % step, on_queue=True)

def plot_all_seeds(prefix):
	import matplotlib.pyplot as plt
	import matplotlib
	#col = [(0.0,0.2,0.8),(1.0,0.1,0.1),(0.0,0.7,0.3),(0.9,0.5,0.1)]
	col = [(0.9047, 0.1918, 0.1988), (0.3718, 0.7176, 0.3612), (0.2941, 0.5447, 0.7494), (1.0000, 0.5482, 0.1000), (0.9550, 0.8946, 0.4722), (0.6859, 0.4035, 0.2412), (0.9718, 0.5553, 0.7741), (0.6400, 0.6400, 0.6400), (0.6365, 0.3753, 0.6753),(0.9047, 0.1918, 0.1988), (0.2941, 0.5447, 0.7494), (0.3718, 0.7176, 0.3612), (1.0000, 0.5482, 0.1000), (0.9550, 0.8946, 0.4722), (0.6859, 0.4035, 0.2412), (0.9718, 0.5553, 0.7741), (0.6400, 0.6400, 0.6400), (0.6365, 0.3753, 0.6753)]
	#good_seeds = [2,1,1,3,1,4,3,1,3]
	good_seeds = [[1,2,3,4,5,6,7,8,9,10],[]]
	font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
	matplotlib.rc('font', **font)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	#ax.set_ylabel("...",labelpad=20)
	for dot2 in range(1,2):
		dot1 = 1
		#for seed in good_seeds[dot2-1]:
		count = 0
		for seed in [1,7]:
			count += 1
		#for seed in good_seeds[dot2-1]:
			#for run in range(1):
			#seed = good_seeds[dot2-1]
			xx = []; yy = []
			for line in open('lammps/'+prefix+'s'+str(seed)+'_'+str(dot1)+'_'+str(dot2)+'.pmf'):
				a = line.split()
				if a:
					try:
						xx.append( float(a[0]) )
						yy.append( float(a[1]) )
					except ValueError: pass
			yy = [y-yy[-1] for y in yy]

			plt.plot( xx, yy, label='run '+str(count),color=col[count+1], linewidth=2.0)
			#plt.gcf().subplots_adjust(left=0.45)
			#plt.plot( xx, yy, marker='.', label=str(seed), color=col[seed-1])
			#leg_text += str(dot2)+'+'+str(dot2)
			#leg_text += str(seed)
			plt.xlabel(r'Radius of Gyration ($\AA$)')
			plt.ylabel('Free Energy (kcal/mol)')
			#plt.legend(loc='lower right')
			#plt.title('metadynamics '+str(dot1)+'+'+str(dot2)+', run '+str(seed))

		plt.legend(loc='lower right')
		plt.show()       	

def calculate_core_volume(seed,dot1,dot2):
	d1 = open('c18_'+str(dot1)+'.xyz')
	N1 = d1.readline()
	N1 = N1.split()
	N1 = int(N1[0])
	d2 = open('c18_'+str(dot1)+'.xyz')
	N2 = d2.readline()
	N2 = N2.split()
	N2 = int(N2[0])
	Ntot = N1 + N2
	f = open('gyr_s'+str(seed)+'_'+str(dot1)+'_'+str(dot2)+'.xyz')
	lines = f.readlines()
	atoms = []
	f.close()
	S_atoms = [[]]
	Pb_atoms = [[]]
	for line in lines:
		parts = line.split()
		if len(parts) == 4:
			atoms += [line]

	for frame in range(len(atoms)/Ntot):
		for line in atoms[frame*Ntot:frame*Ntot+dot1]:
			parts = line.split()	
			x = float(parts[1])
			y = float(parts[2])
			z = float(parts[3])
			S_atoms[frame].append( utils.Struct(x=x, y=y, z=z, element='S') )
		
		for line in atoms[frame*Ntot+dot1:frame*Ntot+2*dot1]:
			parts = line.split()	
			x = float(parts[1])
			y = float(parts[2])
			z = float(parts[3])
			Pb_atoms[frame].append( utils.Struct(x=x, y=y, z=z, element='Pb') )
		
		for line in atoms[frame*Ntot+N1:frame*Ntot+N1+dot2]:
			parts = line.split()	
			x = float(parts[1])
			y = float(parts[2])
			z = float(parts[3])
			S_atoms[frame].append( utils.Struct(x=x, y=y, z=z, element='S') )
		
		for line in atoms[frame*Ntot+N1+dot2:frame*Ntot+N1+2*dot2]:
			parts = line.split()	
			x = float(parts[1])
			y = float(parts[2])
			z = float(parts[3])
			Pb_atoms[frame].append( utils.Struct(x=x, y=y, z=z, element='Pb') )
		
	print S_atoms[0].x
	print Pb_atoms[0].x

	calculate_core_volume(1,1,1)
			
def normalize_radius_of_gyration(seed,dot1,dot2):
	# compute radius of ball with equal volume
	S_rad = 1.80 # angstrom
	Pb_rad = 2.02 # angstrom
	
	# rad
	f = open('gyr_s'+str(seed)+'_'+str(dot1)+'_'+str(dot2)+'.pmf')
	lines = f.readlines()
	f.close()
	

def run_react_dots():
	dot1 = 2
	for dot2 in range(2,3):
		for seed in range(11,12):
			#react_dots(params, dot1, dot2, 'dot'+str(dot1), 'dot'+str(dot2), run_name_prefix='gyr_trunc_s%d_' % seed, random_seed=seed+1, on_queue=True, method='metadynamics')
			react_dots(params, dot1, dot2, 'dot'+str(dot1), 'dot'+str(dot2), run_name_prefix='gyr_trunc_abf_s%d_' % seed, random_seed=seed+1, on_queue=True, method='abf')

run_react_dots()

#prefix = 'gyr_trunc_'
#plot_all_seeds(prefix)



