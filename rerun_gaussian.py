import os, copy, sys, random, cPickle, math
sys.path.append("/fs/home/jms875/Library")
import gaussian, utils, filetypes
import numpy as np

# sim4dot2
def mindot1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_1.log')
		i = 1
		for a in atoms:
			if i >= 23 and i <= 32:
				a.z += -0.5+0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim4out2', atoms)
		else:
			f = open('sim4out2.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot1_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)

	mindot1()

def minligBq1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_1.log')
		i = 1
		for a in atoms:
			if i >= 23 and i <= 32:
				a.z += -0.5+0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq1_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq1()

def mindotBq1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_1.log')
		i = 1
		for a in atoms:
			if i >= 23 and i <= 32:
				a.z += -0.5+0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq1_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq1()

def det1():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_1.log')
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 23 and i <= 32:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet1_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet1_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	det1()

# sim6dot2
def mindot2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 25 and i <= 34:
				a.y += -1.0+0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim6out2', atoms)
		else:
			f = open('sim6out2.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot2_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindot2()

def minligBq2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 25 and i <= 34:
				a.y += -1.0+0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq2_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq2()

def mindotBq2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 25 and i <= 34:
				a.y += -1.0+0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq2_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq2()

def det2():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 25 and i <= 34:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet2_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet2_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	det2()

# sim9dot5
def mindot3():
	for step in range(60):
		mol = utils.Molecule('dot3min.arc')
		theta = 0
		R = np.array([[math.cos(theta), -math.sin(theta), 0.0],[math.sin(theta), math.cos(theta), 0.0],[0.0, 0.0, 1.0]])
		mol.rotate(R)
		atoms = mol.atoms
		i = 1
		for a in atoms:
			if i >= 17 and i <= 26:
				a.z += 0.5-0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim9out5', atoms)
		else:
			f = open('sim9out5.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot3_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)

	mindot3()

def minligBq3():
	for step in range(60):
		mol = utils.Molecule('dot3min.arc')
		theta = 0
		R = np.array([[math.cos(theta), -math.sin(theta), 0.0],[math.sin(theta), math.cos(theta), 0.0],[0.0, 0.0, 1.0]])
		mol.rotate(R)
		atoms = mol.atoms
		i = 1
		for a in atoms:
			if i >= 17 and i <= 26:
				a.z += 0.5-0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq3_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq3()

def mindotBq3():
	for step in range(60):
		mol = utils.Molecule('dot3min.arc')
		theta = 0
		R = np.array([[math.cos(theta), -math.sin(theta), 0.0],[math.sin(theta), math.cos(theta), 0.0],[0.0, 0.0, 1.0]])
		mol.rotate(R)
		atoms = mol.atoms
		i = 1
		for a in atoms:
			if i >= 17 and i <= 26:
				a.z += 0.5-0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq3_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq3()

def det3():
	mol = utils.Molecule('dot3min.arc')
	theta = 0
	R = np.array([[math.cos(theta), -math.sin(theta), 0.0],[math.sin(theta), math.cos(theta), 0.0],[0.0, 0.0, 1.0]])
	mol.rotate(R)
	atoms = mol.atoms
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 17 and i <= 26:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet3_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet3_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	det3()

# sim5dot2
def mindot4():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 49 and i <= 58:
				a.y += -1.0+0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim5out2', atoms)
		else:
			f = open('sim5out2.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot4_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindot4()

def minligBq4():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 49 and i <= 58:
				a.y += -1.0+0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq4_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq4()

def mindotBq4():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 49 and i <= 58:
				a.y += -1.0+0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq4_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq4()

def det4():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 49 and i <= 58:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet4_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet4_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	det4()

# sim10dot
def mindot5():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 131 and i <= 140:
				a.z += -1.0+0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim10out', atoms)
		else:
			f = open('sim10out.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot5_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		
	mindot5()

def minligBq5():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 131 and i <= 140:
				a.z += -1.0+0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq5_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq5()

def mindotBq5():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 131 and i <= 140:
				a.z += -1.0+0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq5_re1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq5()

def det5():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 131 and i <= 140:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet5_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet5_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	det5()

# sim11dot
def mindot6():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 63 and i <= 72:
				a.y += -1.0+0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim11out', atoms)
		else:
			f = open('sim11out.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot6_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)

	mindot6()

def minligBq6():
	for step in range(6,8):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 63 and i <= 72:
				a.y += -1.0+0.1*step
				a.element += '-Bq'
			i += 1

		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq6_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq6()

def mindotBq6():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 63 and i <= 72:
				a.y += -1.0+0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq6_1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq6()

def det6():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 63 and i <= 72:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet6_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet6_1', 'SP SCRF(Solvent=Toluene)', procs=1)
	det6()

# sim10dot2
def mindot5_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 21 and i <= 30:
				a.y += -1.0+0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim10out2', atoms)
		else:
			f = open('sim10out2.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )

		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot5_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)

	mindot5_2()

def minligBq5_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 21 and i <= 30:
				a.y += -1.0+0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq5_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq5_2()

def mindotBq5_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 21 and i <= 30:
				a.y += -1.0+0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq5_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq5_2()

def det5_2():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 21 and i <= 30:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet5_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet5_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	det5_2()

# sim11dot2
def mindot6_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		i = 1
		mag = math.sqrt(6.413663**2+(-6.272733)**2+(3.000239)**2)
		for a in atoms:
			if i >= 73 and i <= 82:
				a.x += -(6.413663)/mag+(6.413663)/mag*0.1*step
				a.y += 6.272733/mag-6.272733/mag*0.1*step
				a.z += -3.000239/mag+3.000239/mag*0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim11out2', atoms)
		else:
			f = open('sim11out2.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot6_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)

	mindot6_2()

def minligBq6_2():
	for step in range(25,26):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		i = 1
		mag = math.sqrt(6.413663**2+(-6.272733)**2+(3.000239)**2)
		for a in atoms:
			if i >= 73 and i <= 82:
				a.x += -(6.413663)/mag+(6.413663)/mag*0.1*step
				a.y += 6.272733/mag-6.272733/mag*0.1*step
				a.z += -3.000239/mag+3.000239/mag*0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq6_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq6_2()

def mindotBq6_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		i = 1
		mag = math.sqrt(6.413663**2+(-6.272733)**2+(3.000239)**2)
		for a in atoms:
			if i >= 73 and i <= 82:
				a.x += -(6.413663)/mag+(6.413663)/mag*0.1*step
				a.y += 6.272733/mag-6.272733/mag*0.1*step
				a.z += -3.000239/mag+3.000239/mag*0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq6_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq6_2()

def det6_2():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 73 and i <= 82:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet6_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet6_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	det6_2()
	
# sim5dot5
def mindot4_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 99 and i <= 108:
				a.y += -1.0+0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim5out5', atoms)
		else:
			f = open('sim5out5.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot4_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dot4_2_'+str(step), 'SP', procs=1)
	mindot4_2()

def minligBq4_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 99 and i <= 108:
				a.y += -1.0+0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq4_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligBq4_2_'+str(step), 'SP', procs=1)
	minligBq4_2()

def mindotBq4_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 99 and i <= 108:
				a.y += -1.0+0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq4_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotBq4_2_'+str(step), 'SP', procs=1)
	mindotBq4_2()

def det4_2():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 99 and i <= 108:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet4_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotdet4_2', 'SP', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet4_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligdet4_2', 'SP', procs=1)
	det4_2()

# sim6dot5
def mindot2_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 35 and i <= 44:
				a.y += 1.0-0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim6out5', atoms)
		else:
			f = open('sim6out5.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
					
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot2_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)

	mindot2_2()

def minligBq2_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 35 and i <= 44:
				a.y += 1.0-0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq2_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq2_2()

def mindotBq2_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 35 and i <= 44:
				a.y += 1.0-0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq2_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq2_2()

def det2_2():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 35 and i <= 44:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet2_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet2_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	det2_2()

# sim9dot2
def mindot3_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_3.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 87 and i <= 96:
				a.z += 0.5-0.1*step
			i += 1
		if step == 0:
			filetypes.write_xyz('sim9out2', atoms)
		else:
			f = open('sim9out2.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot3_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)

	mindot3_2()

def minligBq3_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_3.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 87 and i <= 96:
				a.z += 0.5-0.1*step
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq3_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	minligBq3_2()

def mindotBq3_2():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_3.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i >= 87 and i <= 96:
				a.z += 0.5-0.1*step
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq3_2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
	mindotBq3_2()

def det3_2():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_3.log',check_convergence=False)
	detached_dot = []
	detached_lig = []
	i = 1
	for a in atoms:
		if i >= 87 and i <= 96:
			detached_lig += [a]
		else:
			detached_dot += [a]
		i += 1
	gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet3_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet3_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	det3_2()

# optimize dots in solvent:
def opt_in_sol():
	atoms1 = filetypes.parse_xyz('dot1.xyz')
	gaussian.job(atoms1, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_dot1_opt', 'Opt SCRF(Solvent=Toluene)', procs=1)
	atoms2 = filetypes.parse_xyz('dot2.xyz')
	gaussian.job(atoms2, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_dot2_opt', 'Opt SCRF(Solvent=Toluene)', procs=1)
	atoms3 = filetypes.parse_xyz('dot3.xyz')
	gaussian.job(atoms3, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_dot3_opt', 'Opt SCRF(Solvent=Toluene)', procs=1)
	atoms4 = filetypes.parse_xyz('dot4.xyz')
	gaussian.job(atoms4, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_dot4_opt', 'Opt SCRF(Solvent=Toluene)', procs=1)
	atoms5 = filetypes.parse_xyz('dot5.xyz')
	gaussian.job(atoms5, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_dot5_opt', 'Opt SCRF(Solvent=Toluene)', procs=1)

	opt_in_sol()

def ligBq_in_sol1(filename):
	for ligs in range(5):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (3+10*ligs) and i <= (12+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq11'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	ligBq_in_sol1('dot1_1_opt')
	# observe filename is HSEHsol_11_0... for 1_1 runs
		
def dotBq_in_sol1(filename):
	for ligs in range(5):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (3+10*ligs) or i > (12+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq11'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
#dotBq_in_sol1('dot1_opt')
	dotBq_in_sol1('dot1_1_opt')
	
def ligBq_in_sol2(filename):
	for ligs in range(8):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (5+10*ligs) and i <= (14+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq2'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	ligBq_in_sol2('dot2_opt2')
		
def dotBq_in_sol2(filename):
	for ligs in range(8):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (5+10*ligs) or i > (14+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq2'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	dotBq_in_sol2('dot2_opt2')
	
def ligBq_in_sol3(filename):
	for ligs in range(10):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (7+10*ligs) and i <= (16+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq3'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	ligBq_in_sol3('dot3_opt2')
		
def dotBq_in_sol3(filename):
	for ligs in range(10):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (7+10*ligs) or i > (16+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq3'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	dotBq_in_sol3('dot3_opt2')

def ligBq_in_sol4(filename):
	for ligs in range(12):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (9+10*ligs) and i <= (18+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq4'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	ligBq_in_sol4('dot4_opt2')
		
def dotBq_in_sol4(filename):
	for ligs in range(12):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (9+10*ligs) or i > (18+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq4'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	dotBq_in_sol4('dot4_opt2')
	
def ligBq_in_sol5(filename):
	for ligs in range(14):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (11+10*ligs) and i <= (20+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq5'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	ligBq_in_sol5('dot5_opt3')
		
def dotBq_in_sol5(filename):
	for ligs in range(14):
		energies, newatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		energies, detatoms = gaussian.parse_atoms('gaussian/HSEHsol_'+filename+'.log')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (11+10*ligs) or i > (20+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq5'+'_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	dotBq_in_sol5('dot5_opt3')


def more_opt_in_sol():
	atoms1 = filetypes.parse_xyz('dot1_1.xyz')
	gaussian.job(atoms1, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_dot1_1_opt', 'Opt SCRF(Solvent=Toluene)', procs=1)
	atoms2 = filetypes.parse_xyz('dot1_17.xyz')
	gaussian.job(atoms2, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_dot1_17_opt', 'Opt SCRF(Solvent=Toluene)', procs=1)

	more_opt_in_sol()

def restart_opt_in_sol(dot_size):
	dot = []
	#os.system('cp gaussian/HSEHsol_dot'+str(dot_size)+'_opt.chk gaussian/HSEHsol_dot'+str(dot_size)+'_opt2.chk')
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_dot'+str(dot_size)+'_opt3', 'Geom=AllCheck Opt Guess=Read', procs=1, previous='HSEHsol_dot'+str(dot_size)+'_opt2')
	#gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_dot'+str(dot_size)+'_opt2', 'Geom=AllCheck Opt Guess=Read', procs=1)
	restart_opt_in_sol(4)

# optimized dot runs
def optdot1():
	atoms = filetypes.parse_xyz('dot1_17.xyz')
	gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dot1_min3', 'SP', procs=1)

	optdot1()

def optligBq1():
	for ligs in range(5):
		newatoms = filetypes.parse_xyz('dot1_17.xyz')
		detatoms = filetypes.parse_xyz('dot1_17.xyz')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (3+10*ligs) and i <= (12+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligBq1_min3_'+str(ligs), 'SP', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligdet1_min3_'+str(ligs), 'SP', procs=1)
	optligBq1()
		
def optdotBq1():
	for ligs in range(5):
		newatoms = filetypes.parse_xyz('dot1_17.xyz')
		detatoms = filetypes.parse_xyz('dot1_17.xyz')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (3+10*ligs) or i > (12+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotBq1_min3_'+str(ligs), 'SP', procs=1)
		gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotdet1_min3_'+str(ligs), 'SP', procs=1)
	optdotBq1()

def optdot2():
	atoms = filetypes.parse_xyz('dot2.xyz')
	gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot2_min', 'SP SCRF(Solvent=Toluene)', procs=1)

	optdot2()

def optligBq2():
	for ligs in range(8):
		newatoms = filetypes.parse_xyz('dot2.xyz')
		detatoms = filetypes.parse_xyz('dot2.xyz')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (5+10*ligs) and i <= (14+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq2_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet2_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	optligBq2()
		
def optdotBq2():
	for ligs in range(8):
		newatoms = filetypes.parse_xyz('dot2.xyz')
		detatoms = filetypes.parse_xyz('dot2.xyz')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (5+10*ligs) or i > (14+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq2_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet2_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	optdotBq2()

def optdot3():
	atoms = filetypes.parse_xyz('dot3.xyz')
	gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot3_min', 'SP SCRF(Solvent=Toluene)', procs=1)

	optdot3()

def optligBq3():
	for ligs in range(10):
		newatoms = filetypes.parse_xyz('dot3.xyz')
		detatoms = filetypes.parse_xyz('dot3.xyz')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (7+10*ligs) and i <= (16+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq3_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet3_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	optligBq3()
		
def optdotBq3():
	for ligs in range(10):
		newatoms = filetypes.parse_xyz('dot3.xyz')
		detatoms = filetypes.parse_xyz('dot3.xyz')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (7+10*ligs) or i > (16+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq3_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet3_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	optdotBq3()

def optdot4():
	atoms = filetypes.parse_xyz('dot4.xyz')
	gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot4_min', 'SP SCRF(Solvent=Toluene)', procs=1)

	optdot4()

def optligBq4():
	for ligs in range(12):
		newatoms = filetypes.parse_xyz('dot4.xyz')
		detatoms = filetypes.parse_xyz('dot4.xyz')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (9+10*ligs) and i <= (18+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq4_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet4_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	optligBq4()
		
def optdotBq4():
	for ligs in range(12):
		newatoms = filetypes.parse_xyz('dot4.xyz')
		detatoms = filetypes.parse_xyz('dot4.xyz')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (9+10*ligs) or i > (18+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq4_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet4_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	optdotBq4()

def optdot5():
	atoms = filetypes.parse_xyz('dot5.xyz')
	gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot5_min', 'SP SCRF(Solvent=Toluene)', procs=1)

	optdot5()

def optligBq5():
	for ligs in range(14):
		newatoms = filetypes.parse_xyz('dot5.xyz')
		detatoms = filetypes.parse_xyz('dot5.xyz')
		detached_lig = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i >= (11+10*ligs) and i <= (20+10*ligs):
				detached_lig += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq5_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet5_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	optligBq5()
		
def optdotBq5():
	for ligs in range(14):
		newatoms = filetypes.parse_xyz('dot5.xyz')
		detatoms = filetypes.parse_xyz('dot5.xyz')
		detached_dot = []
		i = 1
		for a,b in zip(newatoms,detatoms):
			if i < (11+10*ligs) or i > (20+10*ligs):
				detached_dot += [b]
				a.element += '-Bq'
			i += 1
		gaussian.job(newatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq5_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet5_min_'+str(ligs), 'SP SCRF(Solvent=Toluene)', procs=1)
	optdotBq5()

# core separation runs
def separate_dot(dot_size,plane):
# note: dot size 2 has only 2 possible planes of separation; dots 3 and 4 have 3 planes of separation
	f = open('/fs/home/ja539/Documents/dot size/minimized dots/dot'+str(dot_size)+'.xyz')
	atoms = f.readlines()
	f.close()
	
	indices1 = []
	indices2 = []
	
	header = 2
	
	dot2_0 = [14,24,44,64]
	dot2_1 = [4,34,54,74]

	dot4_0 = [38,68,118]
	dot4_1 = [8,78,88]
	dot4_2 = [48,58,98]
	dot4_3 = [18,28,108]
	
	if dot_size > 2 and plane == 'xy':
		if dot_size == 3:
			SA = [0,1]
			PbA = [3,5]
			SB = [2]
			PbB = [4]
			dot3_0 = [36,76,96]
			dot3_1 = [6,26,46]
			dot3_2 = [16,56,66,86]
			
			indices1 += SA + PbA
			indices2 += SB + PbB
			
			for a in dot3_0+dot3_1:
				for i in range(10):
					indices1 += [a+i]
			
			for b in dot3_2:
				for i in range(10):
					indices2 += [b+i]
			
		if dot_size == 4:
			SA = [0,3]
			PbA = [4,5]
			SB = [1,2]
			PbB = [6,7]
			
			indices1 += SA + PbA
			indices2 += SB + PbB
				
			for a in dot4_0+dot4_3:
				for i in range(10):
					indices1 += [a+i]
			
			for b in dot4_1+dot4_2:
				for i in range(10):
					indices2 += [b+i]
	
	if plane == 'xz':
		if dot_size == 2:
			SA = [0]
			PbA = [3]
			SB = [1]
			PbB = [2]
			
			indices1 += SA + PbA
			indices2 += SB + PbB
				
			for a in dot2_0:
				for i in range(10):
					indices1 += [a+i]
			
			for b in dot2_1:
				for i in range(10):
					indices2 += [b+i]
			
		if dot_size == 3:
			SA = [1,2]
			PbA = [5,4]
			SB = [0]
			PbB = [3]
			dot3_0 = [36,56,76,96]
			dot3_1 = [6,26,46]
			dot3_2 = [16,66,86]
			
			indices1 += SA + PbA
			indices2 += SB + PbB
			
			for a in dot3_1+dot3_2:
				for i in range(10):
					indices1 += [a+i]
			
			for b in dot3_0:
				for i in range(10):
					indices2 += [b+i]
			
		if dot_size == 4:
			SA = [0,2]
			PbA = [5,6]
			SB = [1,3]
			PbB = [4,7]
			
			indices1 += SA + PbA
			indices2 += SB + PbB
				
			for a in dot4_0+dot4_2:
				for i in range(10):
					indices1 += [a+i]
			
			for b in dot4_1+dot4_3:
				for i in range(10):
					indices2 += [b+i]
	
	if plane == 'yz':
		if dot_size == 2:
			SA = [0]
			PbA = [2]
			SB = [1]
			PbB = [3]
			
			indices1 += SA + PbA
			indices2 += SB + PbB
				
			for a in dot2_0:
				for i in range(10):
					indices1 += [a+i]
			
			for b in dot2_1:
				for i in range(10):
					indices2 += [b+i]
					
		if dot_size == 3:
			SA = [0,2]
			PbA = [3,4]
			SB = [1]
			PbB = [5]
			dot3_0 = [36,76,96]
			dot3_1 = [6,26,66,46]
			dot3_2 = [16,56,86]
			
			indices1 += SA + PbA
			indices2 += SB + PbB
			
			for a in dot3_0+dot3_2:
				for i in range(10):
					indices1 += [a+i]
			
			for b in dot3_1:
				for i in range(10):
					indices2 += [b+i]
	
		if dot_size == 4:
			SA = [0,1]
			PbA = [4,6]
			SB = [2,3]
			PbB = [5,7]
			
			indices1 += SA + PbA
			indices2 += SB + PbB
				
			for a in dot4_0+dot4_1:
				for i in range(10):
					indices1 += [a+i]
			
			for b in dot4_2+dot4_3:
				for i in range(10):
					indices2 += [b+i]
	
	return indices1, indices2
	
def coredot(dot_size, plane):
	indices1, indices2 = separate_dot(dot_size, plane)
	f = open('/fs/home/ja539/Documents/MD_simulations/lammps/spldot_'+str(dot_size)+'_'+plane+'.xyz')
	ref_atoms = f.readlines()
	f.close()
	header = 2
	for step in range(50,51):
		atoms = filetypes.parse_xyz('/fs/home/ja539/Documents/dot size/minimized dots/dot'+str(dot_size)+'.xyz')
		N = len(atoms)
		i = 0
		for a in atoms:
			atom = ref_atoms[header*(step+1) +N*step + i]
			atom = atom.split()
			a.x = float(atom[1])
			a.y = float(atom[2])
			a.z = float(atom[3])
			i += 1
		
		if step == 0:
			filetypes.write_xyz('sim12out'+str(dot_size)+'_'+plane, atoms)
		#else:
			#h = open('sim12out'+str(dot_size)+'_'+plane+'.xyz', 'a')
			#h.write(str(len(atoms))+'\nAtoms\n')
			#for a in atoms:
				#h.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
			
			#h.close()
			
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dot'+str(dot_size)+'_'+str(plane)+'_'+str(step), 'SP', procs=1)
		
	#coredot(2,'yz')
	#coredot(3,'xy')
	coredot(3,'yz')
	#coredot(4,'xy')
	#coredot(4,'yz')

def piece1Bq(dot_size, plane):
	indices1, indices2 = separate_dot(dot_size, plane)
	f = open('/fs/home/ja539/Documents/MD_simulations/lammps/spldot_'+str(dot_size)+'_'+plane+'.xyz')
	ref_atoms = f.readlines()
	f.close()
	header = 2
	for step in range(60):
		atoms = filetypes.parse_xyz('/fs/home/ja539/Documents/dot size/minimized dots/dot'+str(dot_size)+'.xyz')
		N = len(atoms)
		i = 0
		for a in atoms:
			atom = ref_atoms[header*(step+1) +N*step + i]
			#print atom
			atom = atom.split()
			a.x = float(atom[1])
			a.y = float(atom[2])
			a.z = float(atom[3])
			i += 1
		
		for i in indices1:
			atoms[i].element += '-Bq'	
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot'+str(dot_size)+'_'+str(plane)+'_Bq1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		
	piece1Bq(2,'yz')
	piece1Bq(3,'xy')
	piece1Bq(3,'yz')
	piece1Bq(4,'xy')
	piece1Bq(4,'yz')
	
def piece2Bq(dot_size, plane):
	indices1, indices2 = separate_dot(dot_size, plane)
	f = open('/fs/home/ja539/Documents/MD_simulations/lammps/spldot_'+str(dot_size)+'_'+plane+'.xyz')
	ref_atoms = f.readlines()
	f.close()
	header = 2
	for step in range(60):
		atoms = filetypes.parse_xyz('/fs/home/ja539/Documents/dot size/minimized dots/dot'+str(dot_size)+'.xyz')
		N = len(atoms)
		i = 0
		for a in atoms:
			atom = ref_atoms[header*(step+1) +N*step + i]
			#print atom
			atom = atom.split()
			a.x = float(atom[1])
			a.y = float(atom[2])
			a.z = float(atom[3])
			i += 1
		
		for i in indices2:
			atoms[i].element += '-Bq'	
		
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot'+str(dot_size)+'_'+str(plane)+'_Bq2_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		
	piece2Bq(2,'yz')
	piece2Bq(3,'xy')
	piece2Bq(3,'yz')
	piece2Bq(4,'xy')
	piece2Bq(4,'yz')
	
def lone_halves(dot_size, plane):
	indices1, indices2 = separate_dot(dot_size, plane)
	atoms = filetypes.parse_xyz('/fs/home/ja539/Documents/dot size/minimized dots/dot'+str(dot_size)+'.xyz')
	atoms1 = []
	atoms2 = []
	for i in indices1:
		atoms1 += [atoms[i]]
	
	for i in indices2:
		atoms2 += [atoms[i]]
	
	gaussian.job(atoms1, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot'+str(dot_size)+'_piece1_'+str(plane), 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(atoms2, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot'+str(dot_size)+'_piece2_'+str(plane), 'SP SCRF(Solvent=Toluene)', procs=1)
	
	lone_halves(2,'yz')
	lone_halves(3,'xy')
	lone_halves(3,'yz')
	lone_halves(4,'xy')
	lone_halves(4,'yz')
	
# make xyz files

def make_xyz_file(filename, atoms, energy, opt):
	f = open(filename+'.xyz', 'w')
	f.write(str(len(atoms))+'\n')
	f.write(str(energy)+' '+str(opt)+'\n')
	for a in atoms:
		f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
	f.close()
	
# unoptimized runs
def make_xyz_unopt(dot_size,ligN,run,header,count_offset):
	for step in range(5):
		dotBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dotBq'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		ligBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_ligBq'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		#energy = (total_energy - dotBq_energy - ligBq_energy - val)*627.51
		energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
		
		# note: PbS 0
	#		S   1
	#		PbO 2
	#		OO	3
	#		OH  4
	#		HO  5
	#		COO 6
	#		CCO 7
	#		HC  8
	
		opt = 0
		dot_filename = header+'_dot'+str(dot_size)+'_'+str(step+count_offset)
		cutoff = dot_size*2+10*(ligN-1)
		dot_atoms = atoms[:cutoff]+atoms[cutoff+10:]
		for i in range(len(dot_atoms)):
			if dot_atoms[i].element == 'S':
				dot_atoms[i].element = '1'
			elif dot_atoms[i].element == 'Pb':
				if dot_atoms[i+1].element == 'Pb':
					dot_atoms[i].element = '0'
				else:
					dot_atoms[i].element = '2'
			elif dot_atoms[i].element == 'O':
				if dot_atoms[i+1].element == 'H':
					dot_atoms[i].element = '4'
				else:
					dot_atoms[i].element = '3'
			elif dot_atoms[i].element == 'C':
				if dot_atoms[i+1].element == 'H':
					dot_atoms[i].element = '7'
				else:
					dot_atoms[i].element = '6'
			else:
				if dot_atoms[i-1].element == '4':
					dot_atoms[i].element = '5'
				else:
					dot_atoms[i].element = '8'
				
		complex_filename = header+'_complex'+str(dot_size)+'_'+str(step+count_offset)
		complex_atoms = atoms[cutoff:cutoff+10]
		for i in range(len(complex_atoms)):
			if complex_atoms[i].element == 'Pb':
				complex_atoms[i].element = '2'
			elif complex_atoms[i].element == 'O':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '4'
				else:
					complex_atoms[i].element = '3'
			elif complex_atoms[i].element == 'C':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '7'
				else:
					complex_atoms[i].element = '6'
			else:
				if complex_atoms[i-1].element == '4':
					complex_atoms[i].element = '5'
				else:
					complex_atoms[i].element = '8'
		make_xyz_file(dot_filename, dot_atoms, energy, opt)
		make_xyz_file(complex_filename, complex_atoms, energy, opt)

# optimized runs
def make_xyz_opt(dot_size,ligs,run,header,count_offset):
	for step in range(ligs):
		total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'.log')
		dotBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dotBq'+str(dot_size)+'_'+str(step)+'.log')
		ligBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_ligBq'+str(dot_size)+'_'+str(step)+'.log')

		energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
 		opt = 1
		dot_filename = header+'_dot'+str(dot_size)+'_'+str(step+count_offset)
		cutoff = dot_size*2+10*(step)
		dot_atoms = atoms[:cutoff]+atoms[cutoff+10:]
		for i in range(len(dot_atoms)):
			if dot_atoms[i].element == 'S':
				dot_atoms[i].element = '1'
			elif dot_atoms[i].element == 'Pb':
				if dot_atoms[i+1].element == 'Pb':
					dot_atoms[i].element = '0'
				else:
					dot_atoms[i].element = '2'
			elif dot_atoms[i].element == 'O':
				if dot_atoms[i+1].element == 'H':
					dot_atoms[i].element = '4'
				else:
					dot_atoms[i].element = '3'
			elif dot_atoms[i].element == 'C':
				if dot_atoms[i+1].element == 'H':
					dot_atoms[i].element = '7'
				else:
					dot_atoms[i].element = '6'
			else:
				if dot_atoms[i-1].element == '4':
					dot_atoms[i].element = '5'
				else:
					dot_atoms[i].element = '8'
				
		complex_filename = header+'_complex'+str(dot_size)+'_'+str(step+count_offset)
		complex_atoms = atoms[cutoff:cutoff+10]
		for i in range(len(complex_atoms)):
			if complex_atoms[i].element == 'Pb':
				complex_atoms[i].element = '2'
			elif complex_atoms[i].element == 'O':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '4'
				else:
					complex_atoms[i].element = '3'
			elif complex_atoms[i].element == 'C':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '7'
				else:
					complex_atoms[i].element = '6'
			else:
				if complex_atoms[i-1].element == '4':
					complex_atoms[i].element = '5'
				else:
					complex_atoms[i].element = '8'
		make_xyz_file(dot_filename, dot_atoms, energy, opt)
		make_xyz_file(complex_filename, complex_atoms, energy, opt)

	
# core separation runs
def make_xyz_core(dot_size,run,header,count_offset):
	
	for step in range(60):
		dotBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'_Bq1_'+str(step)+'.log')
		ligBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'_Bq2_'+str(step)+'.log')
		total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		#energy = (total_energy - dotBq_energy - ligBq_energy - val)*627.51
		energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
		
		# note: PbS 0
	#		S   1
	#		PbO 2
	#		OO	3
	#		OH  4
	#		HO  5
	#		COO 6
	#		CCO 7
	#		HC  8
		
		indices1, indices2 = separate_dot(dot_size, str(run))
		opt = 0
		dot_filename = header+'_dot'+str(dot_size)+'_'+str(step+count_offset)
		total_energy, atoms1 = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		dot_atoms = [0]*len(indices1)
		for i in range(len(indices1)):
			dot_atoms[i] = atoms1[indices1[i]]

		for i in range(len(dot_atoms)):
			if dot_atoms[i].element == 'S':
				dot_atoms[i].element = '1'
			elif dot_atoms[i].element == 'Pb':
				if dot_atoms[i+1].element == 'Pb':
					dot_atoms[i].element = '0'
				else:
					dot_atoms[i].element = '2'
			elif dot_atoms[i].element == 'O':
				if dot_atoms[i+1].element == 'H':
					dot_atoms[i].element = '4'
				else:
					dot_atoms[i].element = '3'
			elif dot_atoms[i].element == 'C':
				if dot_atoms[i+1].element == 'H':
					dot_atoms[i].element = '7'
				else:
					dot_atoms[i].element = '6'
			else:
				if dot_atoms[i-1].element == '4':
					dot_atoms[i].element = '5'
				else:
					dot_atoms[i].element = '8'
				
		complex_filename = header+'_complex'+str(dot_size)+'_'+str(step+count_offset)
		total_energy, atoms2 = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		complex_atoms = [0]*len(indices2)
		for i in range(len(indices2)):
			complex_atoms[i] = atoms2[indices2[i]]

		for i in range(len(complex_atoms)):
			if complex_atoms[i].element == 'S':
				complex_atoms[i].element = '1'
			elif complex_atoms[i].element == 'Pb':
				if complex_atoms[i+1].element == 'Pb':
					complex_atoms[i].element = '0'
				else:
					complex_atoms[i].element = '2'
			elif complex_atoms[i].element == 'O':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '4'
				else:
					complex_atoms[i].element = '3'
			elif complex_atoms[i].element == 'C':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '7'
				else:
					complex_atoms[i].element = '6'
			else:
				if complex_atoms[i-1].element == '4':
					complex_atoms[i].element = '5'
				else:
					complex_atoms[i].element = '8'
					
		make_xyz_file(dot_filename, dot_atoms, energy, opt)
		make_xyz_file(complex_filename, complex_atoms, energy, opt)

	make_xyz_core(4,'xy','HSEHsol',192)

def make_xyz_grouped(dot_size,run,header,count_offset):
	Ncomplex = [5,8,10,12,14,16]
	
	for step in range(60):
		dotBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dotBq'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		ligBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_ligBq'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		#energy = (total_energy - dotBq_energy - ligBq_energy - val)*627.51
		energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
		
		# note: PbS 0
	#		S   1
	#		PbO 2
	#		OO	3
	#		OH  4
	#		HO  5
	#		COO 6
	#		CCO 7
	#		HC  8
	
		opt = 0
		dot_filename = header+'_dot'+str(dot_size)+'_'+str(step+count_offset)
		dot_atoms = atoms[:dot_size*2]
		for i in range(len(dot_atoms)):
			if dot_atoms[i].element == 'S':
				dot_atoms[i].element = '1'
			else:
				dot_atoms[i].element = '0'
				
		complex_filename = header+'_complex'+str(dot_size)+'_'+str(step+count_offset)
		complex_atoms = atoms[dot_size*2:]
		for i in range(len(complex_atoms)):
			if complex_atoms[i].element == 'Pb':
				complex_atoms[i].element = '2'
			elif complex_atoms[i].element == 'O':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '4'
				else:
					complex_atoms[i].element = '3'
			elif complex_atoms[i].element == 'C':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '7'
				else:
					complex_atoms[i].element = '6'
			else:
				if complex_atoms[i-1].element == '4':
					complex_atoms[i].element = '5'
				else:
					complex_atoms[i].element = '8'
		make_xyz_file(dot_filename, dot_atoms, energy, opt)
		make_xyz_file(complex_filename, complex_atoms, energy, opt)

	make_xyz_grouped(6,'gr1','HSEHsol',120)

def dist(a,b):
	d = ((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2)**(0.5)
	return d
	
def no_BSSE_energies(dot_size,ligN,run,header):
	
	dotdet, atoms = gaussian.parse_atoms('gaussian/'+header+'_dotdet'+str(dot_size)+'_'+str(run)+'.log')
	ligdet, atoms = gaussian.parse_atoms('gaussian/'+header+'_ligdet'+str(dot_size)+'_'+str(run)+'.log')
	E_inf = dotdet+ligdet
	f = open('no_BSSE.txt','a')
	
	for step in range(60):
		total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		energy = (total_energy - E_inf)*627.51
		if step < 59:
			f.write(str(energy)+', ')
		else:
			f.write(str(energy)+'\n')
	
	f.close()
	
	no_BSSE_energies(6,7,2,'HSEH')

def with_BSSE_energies(dot_size,ligN,run,header):
	
	#f = open('with_BSSE.txt','a')
	#g = open('with_BSSE_dist.txt','a')
	
	for step in range(1):
		total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_dot'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		if step == 0:
			b = atoms[dot_size*2+10*(ligN-1)]
			print b.element
			dmin = 1.0e20
			S_near = 0
			for i in range(dot_size):
				a = atoms[i]
				d = dist(a,b)
				print d
				if d < dmin:
					dmin = d
					S_near = i
		print S_near		
		ligBq_energy, atoms1 = gaussian.parse_atoms('gaussian/'+header+'_ligBq'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		dotBq_energy, atoms2 = gaussian.parse_atoms('gaussian/'+header+'_dotBq'+str(dot_size)+'_'+str(run)+'_'+str(step)+'.log')
		d = dist(atoms[S_near],atoms[dot_size*2+10*(ligN-1)])
		energy = (total_energy - ligBq_energy - dotBq_energy)*627.51

	with_BSSE_energies(6,7,2,'HSEH')
		
	
	
# refining parameters runs: selective grouping to make certain interactions internal to the system and limit the parameter count initially, then gradually build parameter set

# group 1: removing all complexes simultaneously
def dot1gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_1.log')
		i = 1
		for a in atoms:
			if i == 3 or i == 13 or i == 23 or i == 33 or i == 43:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 3:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			i += 1
		if step == 0:
			filetypes.write_xyz('dot1gr1', atoms)
		else:
			f = open('dot1gr1.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot1_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dot1_gr1_'+str(step), 'SP', procs=1)
	dot1gr1()
	
def ligBq1gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_1.log')
		i = 1
		for a in atoms:
			if i == 3 or i == 13 or i == 23 or i == 33 or i == 43:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 3:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq1_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligBq1_gr1_'+str(step), 'SP', procs=1)
	ligBq1gr1()
	
def dotBq1gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_1.log')
		i = 1
		for a in atoms:
			if i == 3 or i == 13 or i == 23 or i == 33 or i == 43:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 3:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq1_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotBq1_gr1_'+str(step), 'SP', procs=1)
	dotBq1gr1()

def dot1gr1det():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_1.log')
	dot = atoms[:2]
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet1_gr1_0', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotdet1_gr1_0', 'SP', procs=1)
	for ligs in range(5):
		detatoms = filetypes.parse_xyz('dot1.xyz')
		detached_lig = []
		i = 1
		for b in detatoms:
			if i >= (3+10*ligs) and i <= (12+10*ligs):
				detached_lig += [b]
			i += 1
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet1_gr1_'+str(ligs+1), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligdet1_gr1_'+str(ligs+1), 'SP', procs=1)
	dot1gr1det()
		

def dot2gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 5 or i == 15 or i == 25 or i == 35 or i == 45 or i == 55 or i == 65 or i ==75:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 5:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			i += 1
		if step == 0:
			filetypes.write_xyz('dot2gr1', atoms)
		else:
			f = open('dot2gr1.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot2_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dot2_gr1_'+str(step), 'SP', procs=1)
	dot2gr1()
	
def ligBq2gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 5 or i == 15 or i == 25 or i == 35 or i == 45 or i == 55 or i == 65 or i ==75:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 5:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq2_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligBq2_gr1_'+str(step), 'SP', procs=1)
	ligBq2gr1()
	
def dotBq2gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 5 or i == 15 or i == 25 or i == 35 or i == 45 or i == 55 or i == 65 or i ==75:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 5:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq2_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotBq2_gr1_'+str(step), 'SP', procs=1)
	dotBq2gr1()

def dot2gr1det():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
	dot = atoms[:4]
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet2_gr1_0', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotdet2_gr1_0', 'SP', procs=1)
	for ligs in range(8):
		energies, detatoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_2.log',check_convergence=False)
		detached_lig = []
		i = 1
		for b in detatoms:
			if i >= (5+10*ligs) and i <= (14+10*ligs):
				detached_lig += [b]
			i += 1
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet2_gr1_'+str(ligs+1), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligdet2_gr1_'+str(ligs+1), 'SP', procs=1)
	dot2gr1det()
	
def dot3gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_3.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 7 or i == 17 or i == 27 or i == 37 or i == 47 or i == 57 or i == 67 or i ==77 or i == 87 or i ==97:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 7:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			i += 1
		if step == 0:
			filetypes.write_xyz('dot3gr1', atoms)
		else:
			f = open('dot3gr1.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot3_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dot3_gr1_'+str(step), 'SP', procs=1)
	dot3gr1()
	
def ligBq3gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_3.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 7 or i == 17 or i == 27 or i == 37 or i == 47 or i == 57 or i == 67 or i ==77 or i == 87 or i ==97:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 7:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq3_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligBq3_gr1_'+str(step), 'SP', procs=1)
	ligBq3gr1()
	
def dotBq3gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_3.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 7 or i == 17 or i == 27 or i == 37 or i == 47 or i == 57 or i == 67 or i ==77 or i == 87 or i ==97:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 7:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq3_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotBq3_gr1_'+str(step), 'SP', procs=1)
	dotBq3gr1()

def dot3gr1det():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_3.log',check_convergence=False)
	dot = atoms[:6]
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet3_gr1_0', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotdet3_gr1_0', 'SP', procs=1)
	for ligs in range(10):
		energies, detatoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_3.log',check_convergence=False)
		detached_lig = []
		i = 1
		for b in detatoms:
			if i >= (7+10*ligs) and i <= (16+10*ligs):
				detached_lig += [b]
			i += 1
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet3_gr1_'+str(ligs+1), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligdet3_gr1_'+str(ligs+1), 'SP', procs=1)
	dot3gr1det()

def dot4gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 9 or i == 19 or i == 29 or i == 39 or i == 49 or i == 59 or i == 69 or i ==79 or i == 89 or i ==99 or i == 109 or i == 119:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 9:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			i += 1
		if step == 0:
			filetypes.write_xyz('dot4gr1', atoms)
		else:
			f = open('dot4gr1.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot4_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dot4_gr1_'+str(step), 'SP', procs=1)
	dot4gr1()

def ligBq4gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 9 or i == 19 or i == 29 or i == 39 or i == 49 or i == 59 or i == 69 or i ==79 or i == 89 or i ==99 or i == 109 or i == 119:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 9:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq4_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligBq4_gr1_'+str(step), 'SP', procs=1)
	ligBq4gr1()
	
def dotBq4gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 9 or i == 19 or i == 29 or i == 39 or i == 49 or i == 59 or i == 69 or i ==79 or i == 89 or i ==99 or i == 109 or i == 119:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 9:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq4_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotBq4_gr1_'+str(step), 'SP', procs=1)
	dotBq4gr1()

def dot4gr1det():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
	dot = atoms[:8]
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet4_gr1_0', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotdet4_gr1_0', 'SP', procs=1)
	for ligs in range(12):
		energies, detatoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_4.log',check_convergence=False)
		detached_lig = []
		i = 1
		for b in detatoms:
			if i >= (9+10*ligs) and i <= (18+10*ligs):
				detached_lig += [b]
			i += 1
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet4_gr1_'+str(ligs+1), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligdet4_gr1_'+str(ligs+1), 'SP', procs=1)
	dot4gr1det()
	
def dot5gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 11 or i == 21 or i == 31 or i == 41 or i == 51 or i == 61 or i == 71 or i ==81 or i == 91 or i ==101 or i == 111 or i == 121 or i == 131 or i == 141:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 11:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			i += 1
		if step == 0:
			filetypes.write_xyz('dot5gr1', atoms)
		else:
			f = open('dot5gr1.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot5_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dot5_gr1_'+str(step), 'SP', procs=1)
	dot5gr1()
	
def ligBq5gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 11 or i == 21 or i == 31 or i == 41 or i == 51 or i == 61 or i == 71 or i ==81 or i == 91 or i ==101 or i == 111 or i == 121 or i == 131 or i == 141:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 11:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq5_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligBq5_gr1_'+str(step), 'SP', procs=1)
	ligBq5gr1()
	
def dotBq5gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 11 or i == 21 or i == 31 or i == 41 or i == 51 or i == 61 or i == 71 or i ==81 or i == 91 or i ==101 or i == 111 or i == 121 or i == 131 or i == 141:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 11:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq5_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotBq5_gr1_'+str(step), 'SP', procs=1)
	dotBq5gr1()

def dot5gr1det():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
	dot = atoms[:10]
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet5_gr1_0', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotdet5_gr1_0', 'SP', procs=1)
	for ligs in range(14):
		energies, detatoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_5.log',check_convergence=False)
		detached_lig = []
		i = 1
		for b in detatoms:
			if i >= (11+10*ligs) and i <= (20+10*ligs):
				detached_lig += [b]
			i += 1
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet5_gr1_'+str(ligs+1), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligdet5_gr1_'+str(ligs+1), 'SP', procs=1)
	dot5gr1det()
	
def dot6gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 13 or i == 23 or i == 33 or i == 43 or i == 53 or i == 63 or i == 73 or i ==83 or i == 93 or i ==103 or i == 113 or i == 123 or i == 133 or i == 143 or i == 153 or i == 163:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 13:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			i += 1
		if step == 0:
			filetypes.write_xyz('dot6gr1', atoms)
		else:
			f = open('dot6gr1.xyz','a')
			f.write(str(len(atoms))+'\nAtoms\n')
			for a in atoms:
				f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dot6_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dot6_gr1_'+str(step), 'SP', procs=1)
	dot6gr1()

def ligBq6gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 13 or i == 23 or i == 33 or i == 43 or i == 53 or i == 63 or i == 73 or i ==83 or i == 93 or i ==103 or i == 113 or i == 123 or i == 133 or i == 143 or i == 153 or i == 163:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 13:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligBq6_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligBq6_gr1_'+str(step), 'SP', procs=1)
	ligBq6gr1()
	
def dotBq6gr1():
	for step in range(60):
		energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		i = 1
		for a in atoms:
			if i == 13 or i == 23 or i == 33 or i == 43 or i == 53 or i == 63 or i == 73 or i ==83 or i == 93 or i ==103 or i == 113 or i == 123 or i == 133 or i == 143 or i == 153 or i == 163:
				modulus = (a.x*a.x + a.y*a.y + a.z*a.z)**(0.5)
				xstep = a.x/modulus
				ystep = a.y/modulus
				zstep = a.z/modulus
			if i >= 13:
				a.x += 0.1*(step-5)*xstep
				a.y += 0.1*(step-5)*ystep
				a.z += 0.1*(step-5)*zstep
			else:
				a.element += '-Bq'
			i += 1
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotBq6_gr1_'+str(step), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotBq6_gr1_'+str(step), 'SP', procs=1)
	dotBq6gr1()

def dot6gr1det():
	energies, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
	dot = atoms[:12]
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_dotdet6_gr1_0', 'SP SCRF(Solvent=Toluene)', procs=1)
	gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_dotdet6_gr1_0', 'SP', procs=1)
	for ligs in range(16):
		energies, detatoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		detached_lig = []
		i = 1
		for b in detatoms:
			if i >= (13+10*ligs) and i <= (22+10*ligs):
				detached_lig += [b]
			i += 1
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_ligdet6_gr1_'+str(ligs+1), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(detached_lig, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_ligdet6_gr1_'+str(ligs+1), 'SP', procs=1)
	dot6gr1det()
	
def make_xyz_sol():
	energies, atoms = gaussian.parse_atoms('gaussian/HSEHsol_dot1_1_opt.log')
	filetypes.write_xyz('sol_dot1_1',atoms)

	make_xyz_sol()

def iterative_runs(Niter):
	for i in range(1,5):
		for dot_size in range(1,5):
			# open jovana_test files and take the last geometry
			f = open('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size)+'.xyz')
			Natoms = int(f.readline())
			f.close()
			atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size)+'.xyz')
			atoms = atoms[-Natoms:]
			gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_test_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
			#gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_test_'+str(i)+'_'+str(dot_size), 'SP', procs=1)
				# now isolate the ligands and dots
			Nlig = (Natoms-2*dot_size)/10
			for j in range(Nlig):
				ligatoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size)+'.xyz')
				ligatoms = ligatoms[-Natoms:]
				ligBqatoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size)+'.xyz')
				ligBqatoms = ligBqatoms[-Natoms:]
				dotatoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size)+'.xyz')
				dotatoms = dotatoms[-Natoms:]
				dotBqatoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size)+'.xyz')
				dotBqatoms = dotBqatoms[-Natoms:]
				ligatoms = ligatoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
				ligBqatoms = ligBqatoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
				for a in ligBqatoms:
					a.element += '-Bq'
				dotatoms = dotatoms[0:2*dot_size+j*10]+dotatoms[2*dot_size+(j+1)*10:]
				dotBqatoms = dotBqatoms[0:2*dot_size+j*10]+dotBqatoms[2*dot_size+(j+1)*10:]
				for b in dotBqatoms:
					b.element += '-Bq'
				ligBqatoms += dotatoms
				dotBqatoms += ligatoms
				gaussian.job(ligBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_testligBq_iter'+str(Niter)+'_'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
				gaussian.job(dotBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_testdotBq_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
				#gaussian.job(ligBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_testligBq'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP', procs=1)
				#gaussian.job(dotBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_testdotBq'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP', procs=1)

	iterative_runs(3)

def add_iterative_runs(dot_size,Nruns):
	for i in range(1,Nruns):
		# open jovana_test files and take the last geometry
		f = open('/fs/home/ja539/Desktop/C/jovana_test/jovana_test2_'+str(i)+'_'+str(dot_size)+'.xyz')
		Natoms = int(f.readline())
		f.close()
		atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/jovana_test2_'+str(i)+'_'+str(dot_size)+'.xyz')
		atoms = atoms[-Natoms:]
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_jovana_test2_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
		#gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_test_'+str(i)+'_'+str(dot_size), 'SP', procs=1)
			# now isolate the ligands and dots
		
	add_iterative_runs(1,5)
	
def opt_coord_runs(Niter):
	for i in range(1,11):
		for dot_size in range(4,5):
			# open jovana_test files and take the last geometry
			f = open('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_coord_'+str(i)+'_'+str(dot_size)+'.xyz')
			Natoms = int(f.readline())
			f.close()
			atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_coord_'+str(i)+'_'+str(dot_size)+'.xyz')
			atoms = atoms[-Natoms:]
			gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_coord_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'Opt SCRF(Solvent=Toluene)', procs=1)
			#gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEH_testopt_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'Opt=Loose', procs=1)

	opt_coord_runs(4)

def restart_opt_coord_runs(Niter):
	for i in range(1,11):
		for dot_size in range(1,3):
			dot = []
			gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_coord2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'Geom=AllCheck Opt SCRF(Solvent=Toluene) Guess=Read', procs=1, previous='HSEHsol_coord_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size))

	restart_opt_coord_runs(4)

def opt_iterative_runs(Niter):
	for i in range(1,11):
		for dot_size in range(4,5):
			# open jovana_test files and take the last geometry
			f = open('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size)+'.xyz')
			Natoms = int(f.readline())
			f.close()
			atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size)+'.xyz')
			atoms = atoms[-Natoms:]
			gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_testopt_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'Opt SCRF(Solvent=Toluene)', procs=1)
			#gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEH_testopt_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'Opt=Loose', procs=1)
			
	opt_iterative_runs(4)
		
def restart_opt_iterative_runs(Niter):
	runs = [2, 4, 5, 6, 8]
	for i in runs:
		for dot_size in range(4,5):
			dot = []
			gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEHsol_testopt4_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'Geom=AllCheck Opt SCRF(Solvent=Toluene) Guess=Read', procs=1, previous='HSEHsol_testopt3_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size))
			#gaussian.job(dot, 'HSEH1PBE/LanL2DZ/Auto', 'long', 'HSEH_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'Geom=AllCheck Opt Guess=Read', procs=1, previous='HSEH_testopt_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size))
			#energy, atoms = gaussian.parse_atoms('HSEHsol_jovana_testopt2_'+str(i)+'_'+str(dot_size),check_convergence=False)
			#print energy

	restart_opt_iterative_runs(4)
	
def SP_opt_iterative_runs(Niter):
	runs = [1,3,4]
	for i in runs:
		for dot_size in range(3,4):
			energy, atoms = gaussian.parse_atoms('HSEHsol_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
			#energy, atoms = gaussian.parse_atoms('HSEH_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
			gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_testoptSP'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
			#gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_testoptSPre'+str(Niter)+'_'+str(i)+'_'+str(dot_size), 'SP', procs=1)
				# now isolate the ligands and dots
			Natoms = len(atoms)
			Nlig = (Natoms-2*dot_size)/10
			for j in range(Nlig):
				energy, ligatoms = gaussian.parse_atoms('HSEHsol_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
				energy, ligBqatoms = gaussian.parse_atoms('HSEHsol_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
				energy, dotatoms = gaussian.parse_atoms('HSEHsol_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
				energy, dotBqatoms = gaussian.parse_atoms('HSEHsol_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
				#energy, ligatoms = gaussian.parse_atoms('HSEH_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
				#energy, ligBqatoms = gaussian.parse_atoms('HSEH_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
				#energy, dotatoms = gaussian.parse_atoms('HSEH_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
				#energy, dotBqatoms = gaussian.parse_atoms('HSEH_testopt2_iter'+str(Niter)+'_'+str(i)+'_'+str(dot_size),check_convergence=False)
				ligatoms = ligatoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
				ligBqatoms = ligBqatoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
				for a in ligBqatoms:
					a.element += '-Bq'
				dotatoms = dotatoms[0:2*dot_size+j*10]+dotatoms[2*dot_size+(j+1)*10:]
				dotBqatoms = dotBqatoms[0:2*dot_size+j*10]+dotBqatoms[2*dot_size+(j+1)*10:]
				for b in dotBqatoms:
					b.element += '-Bq'
				ligBqatoms += dotatoms
				dotBqatoms += ligatoms
				gaussian.job(ligBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_testoptSPligBq'+str(Niter)+'_'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
				gaussian.job(dotBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_testoptSPdotBq'+str(Niter)+'_'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
				#gaussian.job(ligBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_testoptSPligBqre'+str(Niter)+'_'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP', procs=1)
				#gaussian.job(dotBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_testoptSPdotBqre'+str(Niter)+'_'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP', procs=1)
				
	SP_opt_iterative_runs(3)
	
def make_xyz_iterative(Niter,count_offset,header='HSEHsol'):
# filein = 'HSEH_jovana_test_'+str(i)+'_'+str(dot_size)
	for dot_size in range(4,5):
		step = 0
		for run in range(1,5):
			# when reading directly from xyz files
			# open jovana_test files and take the last geometry
			#f = open('/fs/home/ja539/Desktop/C/jovana_test/jovana_test_'+str(i)+'_'+str(dot_size)+'.xyz')
			#Natoms = int(f.readline())
			#f.close()
			#atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/jovana_test_'+str(i)+'_'+str(dot_size)+'.xyz')
			#atoms = atoms[-Natoms:]
			#Nlig = (Natoms-2*dot_size)/10
		
			# when reading from log file
			
			energy, atoms = gaussian.parse_atoms('HSEHsol_testopt2_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size),check_convergence=False) # this file only used to obtain Natoms; NOT for the actual energies, as they are continuing to change and no longer compatible with the SP calculations
			Natoms = len(atoms)
			Nlig = (Natoms-2*dot_size)/10
		
			for j in range(Nlig):
				# when reading directly from xyz files
				'''complex_atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/jovana_test_'+str(i)+'_'+str(dot_size)+'.xyz')
				complex_atoms = complex_atoms[-Natoms:]
				dot_atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/jovana_test_'+str(i)+'_'+str(dot_size)+'.xyz')
				dot_atoms = dot_atoms[-Natoms:]
				complex_atoms = complex_atoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
				dot_atoms = dot_atoms[0:2*dot_size+j*10]+dot_atoms[2*dot_size+(j+1)*10:]
			
				dotBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_testdotBq'+str(j)+'_'+str(i)+'_'+str(dot_size)+'.log')
				ligBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_testligBq'+str(j)+'_'+str(i)+'_'+str(dot_size)+'.log')
				total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_test_'+str(i)+'_'+str(dot_size)+'.log')
				energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
			'''
				# when reading from log file
				print j,run,dot_size
				#if ((j==0) and (run==2) and (dot_size==3)) or ((j==5) and (run==3) and (dot_size==3)) or ((j==17) and (run==2) and (dot_size==4)) or ((j==4) and (run==3) and (dot_size==4)) or ((j==18) and (run==4) and (dot_size==4)):
				#if ((j==15) and (run==2) and (dot_size==3)) or ((run==3) and (dot_size==3)) or ((j==14) and (run==4) and (dot_size==3)) or ((run==2) and (dot_size==4)) or ((run==3) and (dot_size==4)) or ((j==3) and (run==4) and (dot_size==2)):
					#step += 1
				if ((run==2) and (dot_size==3)) or ((j==14) and (run==3) and (dot_size==4)) or ((j==15) and (run==3) and (dot_size==4)) or ((j==16) and (run==3) and (dot_size==4)) or ((j==17) and (run==3) and (dot_size==4)) or ((j==18) and (run==3) and (dot_size==4)):
					step = step
				else:
					energy, complex_atoms = gaussian.parse_atoms('HSEHsol_testoptSP'+str(Niter)+'_'+str(run)+'_'+str(dot_size))
					energy, dot_atoms = gaussian.parse_atoms('HSEHsol_testoptSP'+str(Niter)+'_'+str(run)+'_'+str(dot_size))
					complex_atoms = complex_atoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
					dot_atoms = dot_atoms[0:2*dot_size+j*10]+dot_atoms[2*dot_size+(j+1)*10:]
	
					dotBq_energy, atoms = gaussian.parse_atoms('gaussian/HSEHsol_testoptSPdotBq'+str(Niter)+'_'+str(j)+'_'+str(run)+'_'+str(dot_size)+'.log')
					ligBq_energy, atoms = gaussian.parse_atoms('gaussian/HSEHsol_testoptSPligBq'+str(Niter)+'_'+str(j)+'_'+str(run)+'_'+str(dot_size)+'.log')
					total_energy, atoms = gaussian.parse_atoms('gaussian/HSEHsol_testoptSP'+str(Niter)+'_'+str(run)+'_'+str(dot_size)+'.log')
					energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
					
		
			# note: PbS 0 jj
		#		S   1
		#		PbO 2
		#		OO	3
		#		OH  4
		#		HO  5
		#		COO 6
		#		CCO 7
		#		HC  8
	
					opt = 0
					dot_filename = header+'_dot'+str(dot_size)+'_'+str(step+count_offset[dot_size-1])
					for i in range(len(dot_atoms)):
						if dot_atoms[i].element == 'S':
							dot_atoms[i].element = '1'
						elif dot_atoms[i].element == 'Pb':
							if dot_atoms[i+1].element == 'Pb':
								dot_atoms[i].element = '0'
							else:
								dot_atoms[i].element = '2'
						elif dot_atoms[i].element == 'O':
							if dot_atoms[i+1].element == 'H':
								dot_atoms[i].element = '4'
							else:
								dot_atoms[i].element = '3'
						elif dot_atoms[i].element == 'C':
							if dot_atoms[i+1].element == 'H':
								dot_atoms[i].element = '7'
							else:
								dot_atoms[i].element = '6'
						else:
							if dot_atoms[i-1].element == '4':
								dot_atoms[i].element = '5'
							else:
								dot_atoms[i].element = '8'
				
					complex_filename = header+'_complex'+str(dot_size)+'_'+str(step+count_offset[dot_size-1])
					for i in range(len(complex_atoms)):
						if complex_atoms[i].element == 'Pb':
							complex_atoms[i].element = '2'
						elif complex_atoms[i].element == 'O':
							if complex_atoms[i+1].element == 'H':
								complex_atoms[i].element = '4'
							else:
								complex_atoms[i].element = '3'
						elif complex_atoms[i].element == 'C':
							if complex_atoms[i+1].element == 'H':
								complex_atoms[i].element = '7'
							else:
								complex_atoms[i].element = '6'
						else:
							if complex_atoms[i-1].element == '4':
								complex_atoms[i].element = '5'
							else:
								complex_atoms[i].element = '8'
					make_xyz_file(dot_filename, dot_atoms, energy, opt)
					make_xyz_file(complex_filename, complex_atoms, energy, opt)
					step += 1
	count_offset = [223,240,296,252,134]
	make_xyz_iterative(3,count_offset,header='HSEHsol')

def iter_core_runs(dot_size, Niter, run, count_offset):
	#energy, atoms = gaussian.parse_atoms('HSEHsol_testopt2_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size),check_convergence=False)
	#energy, atomsBq = gaussian.parse_atoms('HSEHsol_testopt2_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size),check_convergence=False)
	
	#energy, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
	#energy, atomsBq = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
	Ecore, atoms = gaussian.parse_atoms('HSEHsol_core_dot'+str(dot_size)+'_2')
	Ecore2, atomsBq = gaussian.parse_atoms('HSEHsol_core2_dot'+str(dot_size)+'_2')
	
	coreS = atoms[0:dot_size]
	corePb = atoms[dot_size:2*dot_size]
	S_x = sum([a.z for a in coreS])/dot_size
	Pb_x = sum([a.z for a in corePb])/dot_size
	X = (S_x+Pb_x)/2.0
	core1 = []
	core2 = []
	core1Bq = []
	core2Bq = []
	
	for a in atomsBq:
		a.element += '-Bq'
		
	for i in range(2*dot_size):
		if atoms[i].z > X:
			core1 += [atoms[i]]
			core1Bq += [atomsBq[i]]
		else:
			core2 += [atoms[i]]
			core2Bq += [atomsBq[i]]
	
	for i in range(2*dot_size,len(atoms)):
		if atoms[i].element == 'Pb' and atoms[i].z > X:
			core1 += atoms[i:i+10]
			core1Bq += atomsBq[i:i+10]
		
		if atoms[i].element == 'Pb' and atoms[i].z <= X:
			core2 += atoms[i:i+10]
			core2Bq += atomsBq[i:i+10]
			
	C1 = core1 + core2Bq
	C2 = core2 + core1Bq
	
	#filetypes.write_xyz('HSEHsol_core1_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size), core1)
	#filetypes.write_xyz('HSEHsol_core2_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size), core2)
	filetypes.write_xyz('HSEHsol_core1_dot'+str(dot_size), core1)
	filetypes.write_xyz('HSEHsol_core2_dot'+str(dot_size), core2)
	
	#gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_core_dot'+str(dot_size)+'_3', 'SP SCRF(Solvent=Toluene)', procs=1)
	#gaussian.job(C1, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_core1_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
	#gaussian.job(C2, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_core2_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
	#gaussian.job(C1, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_core1_dot'+str(dot_size)+'_3', 'SP SCRF(Solvent=Toluene)', procs=1)
	#gaussian.job(C2, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_core2_dot'+str(dot_size)+'_3', 'SP SCRF(Solvent=Toluene)', procs=1)
	
	#total_energy, atoms = gaussian.parse_atoms('gaussian/HSEHsol_testoptSP'+str(Niter)+'_'+str(run)+'_'+str(dot_size)+'.log')
	#total_energy, tot_atoms = gaussian.parse_atoms('HSEHsol_core1_dot'+str(dot_size)+'_2')
	
	#Ecore1, atoms = gaussian.parse_atoms('HSEHsol_core1_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size))
	#Ecore2, atoms = gaussian.parse_atoms('HSEHsol_core2_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size))
	
	Ecore1, atoms = gaussian.parse_atoms('HSEHsol_core1_dot'+str(dot_size)+'_2')
	Ecore2, atoms = gaussian.parse_atoms('HSEHsol_core2_dot'+str(dot_size)+'_2')
	Ecore, atoms = gaussian.parse_atoms('HSEHsol_core_dot'+str(dot_size)+'_2')
	#for a,b in zip(atoms,tot_atoms):
		#b.x = a.x
		#b.y = a.y
		#b.z = a.z
	#filetypes.write_xyz('out',tot_atoms)
	#gaussian.job(tot_atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_core_dot'+str(dot_size)+'_2', 'SP SCRF(Solvent=Toluene)', procs=1)
	#diff = 0.0
	#for a,b in zip(tot_atoms,atoms):
	#	print a.x, b.x
	#print total_energy, Ecore1, Ecore2
	energy = (Ecore - Ecore1 - Ecore2)*627.51
	dot_atoms = core1
	complex_atoms = core2
	header = 'HSEHsol'
	opt = 0
	dot_filename = header+'_dot'+str(dot_size)+'_'+str(count_offset)
	for i in range(len(dot_atoms)):
		if dot_atoms[i].element == 'S':
			dot_atoms[i].element = '1'
		elif dot_atoms[i].element == 'Pb':
			if dot_atoms[i+1].element == 'Pb':
				dot_atoms[i].element = '0'
			else:
				dot_atoms[i].element = '2'
		elif dot_atoms[i].element == 'O':
			if dot_atoms[i+1].element == 'H':
				dot_atoms[i].element = '4'
			else:
				dot_atoms[i].element = '3'
		elif dot_atoms[i].element == 'C':
			if dot_atoms[i+1].element == 'H':
				dot_atoms[i].element = '7'
			else:
				dot_atoms[i].element = '6'
		else:
			if dot_atoms[i-1].element == '4':
				dot_atoms[i].element = '5'
			else:
				dot_atoms[i].element = '8'
			
	complex_filename = header+'_complex'+str(dot_size)+'_'+str(count_offset)
	for i in range(len(complex_atoms)):
		if complex_atoms[i].element == 'S':
			complex_atoms[i].element = '1'
		elif complex_atoms[i].element == 'Pb':
			if complex_atoms[i+1].element == 'Pb':
				complex_atoms[i].element = '0'
			else:
				complex_atoms[i].element = '2'
		elif complex_atoms[i].element == 'O':
			if complex_atoms[i+1].element == 'H':
				complex_atoms[i].element = '4'
			else:
				complex_atoms[i].element = '3'
		elif complex_atoms[i].element == 'C':
			if complex_atoms[i+1].element == 'H':
				complex_atoms[i].element = '7'
			else:
				complex_atoms[i].element = '6'
		else:
			if complex_atoms[i-1].element == '4':
				complex_atoms[i].element = '5'
			else:
				complex_atoms[i].element = '8'
	#make_xyz_file(dot_filename, dot_atoms, energy, opt)
	#make_xyz_file(complex_filename, complex_atoms, energy, opt)
	#iter_core_runs(4,3,1,323)
	#iter_core_runs(4,3,2,324)
	iter_core_runs(6,3,4,120)

def iter_core_runs2(dot_size, count_offset):
	for step in range(0,1):
		print step
	#energy, atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
	#energy, atomsBq = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		Ecore1, complex_atoms = gaussian.parse_atoms('gaussian/HSEHsol_core1_dot'+str(dot_size)+'_4_'+str(step)+'.log')
		Ecore2, dot_atoms = gaussian.parse_atoms('gaussian/HSEHsol_core2_dot'+str(dot_size)+'_4_'+str(step)+'.log')
		Ecore, atoms = gaussian.parse_atoms('gaussian/HSEHsol_core_dot'+str(dot_size)+'_4_'+str(step)+'.log')
		#Ecore2, atomsBq = gaussian.parse_atoms('gaussian/HSEHsol_core_dot'+str(dot_size)+'_4_'+str(step)+'.log')
	
		energy = (Ecore - Ecore1 - Ecore2)*627.51
		complex_atoms = complex_atoms[0:96]
		dot_atoms = dot_atoms[0:76]
		#print len(core1)
		#print len(core2)
		header = 'HSEHsol'
		opt = 0
		dot_filename = header+'_dot'+str(dot_size)+'_'+str(count_offset+step)
		for i in range(len(dot_atoms)):
			if dot_atoms[i].element == 'S':
				dot_atoms[i].element = '1'
			elif dot_atoms[i].element == 'Pb':
				if dot_atoms[i+1].element == 'Pb':
					dot_atoms[i].element = '0'
				else:
					dot_atoms[i].element = '2'
			elif dot_atoms[i].element == 'O':
				if dot_atoms[i+1].element == 'H':
					dot_atoms[i].element = '4'
				else:
					dot_atoms[i].element = '3'
			elif dot_atoms[i].element == 'C':
				if dot_atoms[i+1].element == 'H':
					dot_atoms[i].element = '7'
				else:
					dot_atoms[i].element = '6'
			else:
				if dot_atoms[i-1].element == '4':
					dot_atoms[i].element = '5'
				else:
					dot_atoms[i].element = '8'
			
		complex_filename = header+'_complex'+str(dot_size)+'_'+str(step+count_offset)
		for i in range(len(complex_atoms)):
			if complex_atoms[i].element == 'S':
				complex_atoms[i].element = '1'
			elif complex_atoms[i].element == 'Pb':
				if complex_atoms[i+1].element == 'Pb':
					complex_atoms[i].element = '0'
				else:
					complex_atoms[i].element = '2'
			elif complex_atoms[i].element == 'O':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '4'
				else:
					complex_atoms[i].element = '3'
			elif complex_atoms[i].element == 'C':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '7'
				else:
					complex_atoms[i].element = '6'
			else:
				if complex_atoms[i-1].element == '4':
					complex_atoms[i].element = '5'
				else:
					complex_atoms[i].element = '8'
		make_xyz_file(dot_filename, dot_atoms, energy, opt)
		make_xyz_file(complex_filename, complex_atoms, energy, opt)
	#iter_core_runs2(4,3,1,323)
	#iter_core_runs2(4,3,2,324)
	iter_core_runs2(6, 133)


def get_absolute_energies(dot_size,Niter,header='HSEHsol'):
	E = []
	min_energy = 1e20
	for dot_size in range(1,2):
		for run in range(1,5):
			energy, atoms = gaussian.parse_atoms('jovana_testoptSP'+str(Niter)+'_'+str(run)+'_'+str(dot_size),check_convergence=False)
			E += [energy]
			if energy < min_energy:
				min_energy = energy
	for i in range(1,11):
		opt_energy, atoms = gaussian.parse_atoms('HSEHsol_jovana_testopt2_'+str(i)+'_'+str(dot_size),check_convergence=False)
		E += [opt_energy]
		if opt_energy < min_energy:
			min_energy = opt_energy
	
	E = [(e-min_energy)*627.5 for e in E]
	print e
	p = [0.0]*len(E)
	for i in range(len(E)):
		p[i] = math.exp(-E[i]/0.8)
		
	import matplotlib.pyplot as plt
	plt.plot(E)
	plt.show()

	get_absolute_energies(1,2)

def make_HSEH_xyz_iterative(Niter,count_offset,header='HSEH'):
# filein = 'HSEH_jovana_test_'+str(i)+'_'+str(dot_size)
	for dot_size in range(3,5):
		step = 0
		for run in range(2,6):
			# when reading directly from xyz files
			# open jovana_test files and take the last geometry
			'''f = open('/fs/home/ja539/Desktop/C/jovana_test/jovana_test_'+str(i)+'_'+str(dot_size)+'.xyz')
			Natoms = int(f.readline())
			f.close()
			atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/jovana_test_'+str(i)+'_'+str(dot_size)+'.xyz')
			atoms = atoms[-Natoms:]
			Nlig = (Natoms-2*dot_size)/10
		'''
			# when reading from log file
			
			energy, atoms = gaussian.parse_atoms('HSEH_testopt2_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size),check_convergence=False) # this file only used to obtain Natoms; NOT for the actual energies, as they are continuing to change and no longer compatible with the SP calculations
			Natoms = len(atoms)
			Nlig = (Natoms-2*dot_size)/10
		
			for j in range(Nlig):
				# when reading directly from xyz files
				'''complex_atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/jovana_test_'+str(i)+'_'+str(dot_size)+'.xyz')
				complex_atoms = complex_atoms[-Natoms:]
				dot_atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/jovana_test_'+str(i)+'_'+str(dot_size)+'.xyz')
				dot_atoms = dot_atoms[-Natoms:]
				complex_atoms = complex_atoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
				dot_atoms = dot_atoms[0:2*dot_size+j*10]+dot_atoms[2*dot_size+(j+1)*10:]
			
				dotBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_testdotBq'+str(j)+'_'+str(i)+'_'+str(dot_size)+'.log')
				ligBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_testligBq'+str(j)+'_'+str(i)+'_'+str(dot_size)+'.log')
				total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_test_'+str(i)+'_'+str(dot_size)+'.log')
				energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
			'''
				# when reading from log file
				print j,run,dot_size
				#if ((j==0) and (run==2) and (dot_size==3)) or ((j==5) and (run==3) and (dot_size==3)) or ((j==17) and (run==2) and (dot_size==4)) or ((j==4) and (run==3) and (dot_size==4)) or ((j==18) and (run==4) and (dot_size==4)):
				#if ((j==0) and (run==2) and (dot_size==1)) or ((run==2) and (dot_size==2)) or ((j==1) and (run==1) and (dot_size==4)) or ((run==2) and (dot_size==4)) or ((run==3) and (dot_size==4)) or ((run==4) and (dot_size==4))  or ((j==18) and (run==5) and (dot_size==4)):
				if ((j==7) and (run==2) and (dot_size==3)) or ((run==4) and (dot_size==3)) or ((j==1) and (run==5) and (dot_size==3)) or ((run==2) and (dot_size==4)):
					step = step
				else:
					energy, complex_atoms = gaussian.parse_atoms('HSEH_testoptSPre'+str(Niter)+'_'+str(run)+'_'+str(dot_size))
					#energy, complex_atoms = gaussian.parse_atoms('HSEH_testoptSP_'+str(run)+'_'+str(dot_size))
					energy, dot_atoms = gaussian.parse_atoms('HSEH_testoptSPre'+str(Niter)+'_'+str(run)+'_'+str(dot_size))
					#energy, dot_atoms = gaussian.parse_atoms('HSEH_testoptSP_'+str(run)+'_'+str(dot_size))
					complex_atoms = complex_atoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
					dot_atoms = dot_atoms[0:2*dot_size+j*10]+dot_atoms[2*dot_size+(j+1)*10:]
	
					dotBq_energy, atoms = gaussian.parse_atoms('gaussian/HSEH_testoptSPdotBqre'+str(Niter)+'_'+str(j)+'_'+str(run)+'_'+str(dot_size)+'.log')
					ligBq_energy, atoms = gaussian.parse_atoms('gaussian/HSEH_testoptSPligBqre'+str(Niter)+'_'+str(j)+'_'+str(run)+'_'+str(dot_size)+'.log')
					total_energy, atoms = gaussian.parse_atoms('gaussian/HSEH_testoptSPre'+str(Niter)+'_'+str(run)+'_'+str(dot_size)+'.log')
					#total_energy, atoms = gaussian.parse_atoms('gaussian/HSEH_testoptSP_'+str(run)+'_'+str(dot_size)+'.log')
					energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
					
		
			# note: PbS 0
		#		S   1
		#		PbO 2
		#		OO	3
		#		OH  4
		#		HO  5
		#		COO 6
		#		CCO 7
		#		HC  8
	
					opt = 0
					dot_filename = header+'_dot'+str(dot_size)+'_'+str(step+count_offset[dot_size-1])
					for i in range(len(dot_atoms)):
						if dot_atoms[i].element == 'S':
							dot_atoms[i].element = '1'
						elif dot_atoms[i].element == 'Pb':
							if dot_atoms[i+1].element == 'Pb':
								dot_atoms[i].element = '0'
							else:
								dot_atoms[i].element = '2'
						elif dot_atoms[i].element == 'O':
							if dot_atoms[i+1].element == 'H':
								dot_atoms[i].element = '4'
							else:
								dot_atoms[i].element = '3'
						elif dot_atoms[i].element == 'C':
							if dot_atoms[i+1].element == 'H':
								dot_atoms[i].element = '7'
							else:
								dot_atoms[i].element = '6'
						else:
							if dot_atoms[i-1].element == '4':
								dot_atoms[i].element = '5'
							else:
								dot_atoms[i].element = '8'
				
					complex_filename = header+'_complex'+str(dot_size)+'_'+str(step+count_offset[dot_size-1])
					for i in range(len(complex_atoms)):
						if complex_atoms[i].element == 'Pb':
							complex_atoms[i].element = '2'
						elif complex_atoms[i].element == 'O':
							if complex_atoms[i+1].element == 'H':
								complex_atoms[i].element = '4'
							else:
								complex_atoms[i].element = '3'
						elif complex_atoms[i].element == 'C':
							if complex_atoms[i+1].element == 'H':
								complex_atoms[i].element = '7'
							else:
								complex_atoms[i].element = '6'
						else:
							if complex_atoms[i-1].element == '4':
								complex_atoms[i].element = '5'
							else:
								complex_atoms[i].element = '8'
					make_xyz_file(dot_filename, dot_atoms, energy, opt)
					make_xyz_file(complex_filename, complex_atoms, energy, opt)
					step += 1
	count_offset = [160,188,250,252,134]
	make_HSEH_xyz_iterative(3,count_offset,header='HSEH')

def parse_xyz_iterative(dot_size,header):
	energy = []
	oenergy = []
	min_energy = 1.0e20
	for i in range(1,11):
		total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_test_'+str(i)+'_'+str(dot_size)+'.log')
		opt_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_testopt_'+str(i)+'_'+str(dot_size)+'.log',check_convergence=False)
		if total_energy < min_energy:
			min_energy = total_energy
		if opt_energy < min_energy:
			min_energy = opt_energy
		
		energy += [total_energy]
		oenergy += [opt_energy]
	
	energy = [(e-min_energy)*627.5 for e in energy]
	oenergy = [(e-min_energy)*627.5 for e in oenergy]
	
	import matplotlib.pyplot as plt
	plt.plot(energy)
	plt.plot(oenergy)
	plt.show()
	
	parse_xyz_iterative(1,'HSEHsol')

def param_validation_runs(Nruns):
	for dot_size in range(1,7):
		for i in range(Nruns):
			# open jovana_test files and take the last geometry
			f = open('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(i)+'_'+str(dot_size)+'.xyz')
			Natoms = int(f.readline())
			f.close()
			atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(i)+'_'+str(dot_size)+'.xyz')
			atoms = atoms[-Natoms:]
			#gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_test_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
			#gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_test_'+str(i)+'_'+str(dot_size), 'SP', procs=1)
				# now isolate the ligands and dots
			Nlig = (Natoms-2*dot_size)/10
			print Nlig
			for j in range(Nlig):
				ligatoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(i)+'_'+str(dot_size)+'.xyz')
				ligatoms = ligatoms[-Natoms:]
				ligBqatoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(i)+'_'+str(dot_size)+'.xyz')
				ligBqatoms = ligBqatoms[-Natoms:]
				dotatoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(i)+'_'+str(dot_size)+'.xyz')
				dotatoms = dotatoms[-Natoms:]
				dotBqatoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(i)+'_'+str(dot_size)+'.xyz')
				dotBqatoms = dotBqatoms[-Natoms:]
				ligatoms = ligatoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
				ligBqatoms = ligBqatoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
				for a in ligBqatoms:
					a.element += '-Bq'
				dotatoms = dotatoms[0:2*dot_size+j*10]+dotatoms[2*dot_size+(j+1)*10:]
				dotBqatoms = dotBqatoms[0:2*dot_size+j*10]+dotBqatoms[2*dot_size+(j+1)*10:]
				for b in dotBqatoms:
					b.element += '-Bq'
				ligBqatoms += dotatoms
				dotBqatoms += ligatoms
				#gaussian.job(ligBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_testligBqre'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
				#gaussian.job(dotBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_testdotBq'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP SCRF(Solvent=Toluene)', procs=1)
				#gaussian.job(ligBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_testligBq'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP', procs=1)
				#gaussian.job(dotBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_testdotBq'+str(j)+'_'+str(i)+'_'+str(dot_size), 'SP', procs=1)
				
	param_validation_runs(5)

	
def make_xyz_validation(Nruns,count_offset,header='HSEHsol'):
	for dot_size in range(4):
		step = 0
		for run in range(Nruns):
			# open jovana_test files and take the last geometry
			f = open('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(run)+'_'+str(dot_size)+'.xyz')
			Natoms = int(f.readline())
			f.close()
			atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(run)+'_'+str(dot_size)+'.xyz')
			atoms = atoms[-Natoms:]
			Nlig = (Natoms-2*dot_size)/10
			for j in range(Nlig):
				complex_atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(run)+'_'+str(dot_size)+'.xyz')
				complex_atoms = complex_atoms[-Natoms:]
				dot_atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_'+str(run)+'_'+str(dot_size)+'.xyz')
				dot_atoms = dot_atoms[-Natoms:]
				complex_atoms = complex_atoms[2*dot_size+j*10:2*dot_size+(j+1)*10]
				dot_atoms = dot_atoms[0:2*dot_size+j*10]+dot_atoms[2*dot_size+(j+1)*10:]
				if not (run == 2 and dot_size == 1) and not (j == 2 and run == 3 and dot_size == 1) and not (run == 0 and dot_size == 2) and not (j == 3 and run == 3 and dot_size == 4) and not (j == 7 and run == 2 and dot_size == 2) and not (run == 1 and dot_size == 3) and not (dot_size == 5):
					print j,run,dot_size
					dotBq_energy, atoms = gaussian.parse_atoms('gaussian/HSEHsol_testdotBq'+str(j)+'_'+str(run)+'_'+str(dot_size)+'.log')
					ligBq_energy, atoms = gaussian.parse_atoms('gaussian/HSEHsol_testligBq'+str(j)+'_'+str(run)+'_'+str(dot_size)+'.log')
					total_energy, atoms = gaussian.parse_atoms('gaussian/HSEHsol_test_'+str(run)+'_'+str(dot_size)+'.log')
					energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
		
				# note: PbS 0
			#		S   1
			#		PbO 2
			#		OO	3
			#		OH  4
			#		HO  5
			#		COO 6
			#		CCO 7
			#		HC  8
	
					opt = 0
					dot_filename = header+'_dot'+str(dot_size)+'_'+str(step+count_offset[dot_size-1])
					for i in range(len(dot_atoms)):
						if dot_atoms[i].element == 'S':
							dot_atoms[i].element = '1'
						elif dot_atoms[i].element == 'Pb':
							if dot_atoms[i+1].element == 'Pb':
								dot_atoms[i].element = '0'
							else:
								dot_atoms[i].element = '2'
						elif dot_atoms[i].element == 'O':
							if dot_atoms[i+1].element == 'H':
								dot_atoms[i].element = '4'
							else:
								dot_atoms[i].element = '3'
						elif dot_atoms[i].element == 'C':
							if dot_atoms[i+1].element == 'H':
								dot_atoms[i].element = '7'
							else:
								dot_atoms[i].element = '6'
						else:
							if dot_atoms[i-1].element == '4':
								dot_atoms[i].element = '5'
							else:
								dot_atoms[i].element = '8'
				
					complex_filename = header+'_complex'+str(dot_size)+'_'+str(step+count_offset[dot_size-1])
					for i in range(len(complex_atoms)):
						if complex_atoms[i].element == 'Pb':
							complex_atoms[i].element = '2'
						elif complex_atoms[i].element == 'O':
							if complex_atoms[i+1].element == 'H':
								complex_atoms[i].element = '4'
							else:
								complex_atoms[i].element = '3'
						elif complex_atoms[i].element == 'C':
							if complex_atoms[i+1].element == 'H':
								complex_atoms[i].element = '7'
							else:
								complex_atoms[i].element = '6'
						else:
							if complex_atoms[i-1].element == '4':
								complex_atoms[i].element = '5'
							else:
								complex_atoms[i].element = '8'
					make_xyz_file(dot_filename, dot_atoms, energy, opt)
					make_xyz_file(complex_filename, complex_atoms, energy, opt)
					step += 1
	
	count_offset = [0,0,0,0,0,0] # will go into separate HSEHsol_validation folder; these first ones are validating the original HSEHsol_test3_2_out5 params
	make_xyz_validation(5,count_offset,header='HSEHsol')

def all_grouped_SPOCs():
	for i in range(1,6):
		atoms = filetypes.parse_xyz('dot'+str(i)+'.xyz')
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_jovana_grp_'+str(i), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(atoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_grp_'+str(i), 'SP', procs=1)
		# now isolate the ligands and dots
		Ncore = i*2
		ligatoms = filetypes.parse_xyz('dot'+str(i)+'.xyz')
		ligatoms = ligatoms[Ncore:]
		ligBqatoms = filetypes.parse_xyz('dot'+str(i)+'.xyz')
		ligBqatoms = ligBqatoms[Ncore:]
		dotatoms = filetypes.parse_xyz('dot'+str(i)+'.xyz')
		dotatoms = dotatoms[0:Ncore]
		dotBqatoms = filetypes.parse_xyz('dot'+str(i)+'.xyz')
		dotBqatoms = dotBqatoms[0:Ncore]
		for a in ligBqatoms:
			a.element += '-Bq'
		for b in dotBqatoms:
			b.element += '-Bq'
		ligBqatoms += dotatoms
		dotBqatoms += ligatoms
		gaussian.job(ligBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_jovana_grpligBq'+str(i), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(dotBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEHsol_jovana_grpdotBq'+str(i), 'SP SCRF(Solvent=Toluene)', procs=1)
		gaussian.job(ligBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_grpligBq'+str(i), 'SP', procs=1)
		gaussian.job(dotBqatoms, 'HSEH1PBE/LanL2DZ/Auto', 'batch', 'HSEH_jovana_grpdotBq'+str(i), 'SP', procs=1)

	all_grouped_SPOCs()

def make_xyz_grouped(count_offset,header='HSEHsol'):
# count offset needs to be a vector that's unique for each dot size
# filein = 'HSEH_jovana_grp_'+str(i)
	for dot_size in range(1,6): # all dot sizes
		# open jovana_test_grp files
		Ncore = dot_size*2
		atoms = filetypes.parse_xyz('dot'+str(dot_size)+'.xyz')
		complex_atoms = filetypes.parse_xyz('dot'+str(dot_size)+'.xyz')
		complex_atoms = complex_atoms[Ncore:]
		dot_atoms = filetypes.parse_xyz('dot'+str(dot_size)+'.xyz')
		dot_atoms = dot_atoms[0:Ncore]
	
		dotBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_grpdotBq'+str(dot_size)+'.log')
		ligBq_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_grpligBq'+str(dot_size)+'.log')
		total_energy, atoms = gaussian.parse_atoms('gaussian/'+header+'_jovana_grp_'+str(dot_size)+'.log')
		#energy = (total_energy - dotBq_energy - ligBq_energy - val)*627.51
		energy = (total_energy - dotBq_energy - ligBq_energy)*627.51
		
		# note: PbS 0
	#		S   1
	#		PbO 2
	#		OO	3
	#		OH  4
	#		HO  5
	#		COO 6
	#		CCO 7
	#		HC  8
	
		opt = 1
		dot_filename = header+'_dot'+str(dot_size)+'_'+str(count_offset[dot_size-1])
		for i in range(len(dot_atoms)):
			if dot_atoms[i].element == 'S':
				dot_atoms[i].element = '1'
			else:
				dot_atoms[i].element == '0'
				
		complex_filename = header+'_complex'+str(dot_size)+'_'+str(count_offset[dot_size-1])
		for i in range(len(complex_atoms)):
			if complex_atoms[i].element == 'Pb':
				complex_atoms[i].element = '2'
			elif complex_atoms[i].element == 'O':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '4'
				else:
					complex_atoms[i].element = '3'
			elif complex_atoms[i].element == 'C':
				if complex_atoms[i+1].element == 'H':
					complex_atoms[i].element = '7'
				else:
					complex_atoms[i].element = '6'
			else:
				if complex_atoms[i-1].element == '4':
					complex_atoms[i].element = '5'
				else:
					complex_atoms[i].element = '8'
		make_xyz_file(dot_filename, dot_atoms, energy, opt)
		make_xyz_file(complex_filename, complex_atoms, energy, opt)
		
	count_offset = [70, 188, 250, 252, 134]
	make_xyz_grouped(count_offset,header='HSEHsol')
