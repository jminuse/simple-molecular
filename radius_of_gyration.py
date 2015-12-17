import os, copy, sys, random, cPickle, math
sys.path.append("/fs/home/jms875/Library")
import utils, filetypes, lammps

def center_of_geometry(atoms):
	return utils.Struct(  x = sum(a.x for a in atoms)/len(atoms),
		y = sum(a.y for a in atoms)/len(atoms),
		z = sum(a.z for a in atoms)/len(atoms)   )

def radius_of_gyration(atoms):
		center = center_of_geometry(atoms)
		return math.sqrt( 1.0/len(atoms) * sum( [utils.dist_squared(a,center) for a in atoms] ) )

for dot2 in range(1,2):
	filename = 'gyr_trunc_s7_1_'+str(dot2)
	timesteps = 4008
	#print timesteps
	core_atoms = [[] for _ in xrange(timesteps)]
	t = -1
	count = 0
	f = open(filename+'.xyz')
	g = open(filename+'_rad.xyz','w')
	contents = f.readlines()
	f.close()
	for line in contents: # isolate only core atoms
		line = line.split()
		if line[0] == 'Atoms.':
			t += 1
		if line[0] == 'S':
			#print t, len(line)
			count += 1
			core_atoms[t].append( utils.Struct(x=float(line[1]), y=float(line[2]), z=float(line[3]), element='S') )
		while line[0] != 'S' and (count > 0):
			core_atoms[t].append( utils.Struct(x=float(line[1]), y=float(line[2]), z=float(line[3]), element='Pb') )
			count -= 1
	for t in range(timesteps):
		if len(core_atoms[t]) > 0:
			R = radius_of_gyration(core_atoms[t])
			g.write(str(R)+'\n')
	
	g.close()
