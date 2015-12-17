import os, copy, sys, random, cPickle, math
sys.path.append("/fs/home/jms875/Library")
import gaussian, utils, filetypes
import numpy as np

def make_xyz_file(filename, atoms, action):
	f = open(filename+'.xyz', action)
	f.write(str(len(atoms))+'\nAtoms\n')
	for a in atoms:
		f.write('%s %f %f %f\n' % (a.element, a.x, a.y, a.z) )
	f.close()
	
def weighted_COM(start_atoms,end_atoms):
	MASS_BY_TYPE = {'Pb': 207.2, 'S': 32.065, 'O': 15.9994, 'H': 1.0079, 'C': 12.0107}
# shift all structures so that origin is the COM:
	start_comx = 0.0
	start_comy = 0.0
	start_comz = 0.0
	start_mtot = 0.0
	end_comx = 0.0
	end_comy = 0.0
	end_comz = 0.0
	end_mtot = 0.0
	# calculate center of mass
	for s,e in zip(start_atoms, end_atoms):
		start_comx  += s.x*MASS_BY_TYPE[ s.element ]
		start_comy  += s.y*MASS_BY_TYPE[ s.element ]
		start_comz  += s.z*MASS_BY_TYPE[ s.element ]
		start_mtot += MASS_BY_TYPE[ s.element ]
		
		end_comx  += e.x*MASS_BY_TYPE[ e.element ]
		end_comy  += e.y*MASS_BY_TYPE[ e.element ]
		end_comz  += e.z*MASS_BY_TYPE[ e.element ]
		end_mtot += MASS_BY_TYPE[ e.element ]
		
		start_comx = start_comx/start_mtot
		start_comy = start_comy/start_mtot
		start_comz = start_comz/start_mtot
		
		end_comx = end_comx/end_mtot
		end_comy = end_comy/end_mtot
		end_comz = end_comz/end_mtot

	startx = [a.x-start_comx for a in start_atoms]
	starty = [a.y-start_comy for a in start_atoms]
	startz = [a.z-start_comz for a in start_atoms]

	endx = [a.x-end_comx for a in end_atoms]
	endy = [a.y-end_comy for a in end_atoms]
	endz = [a.z-end_comz for a in end_atoms]
	
		
	for i in range(len(start_atoms)):
		start_atoms[i].x = startx[i]
		start_atoms[i].y = starty[i]
		start_atoms[i].z = startz[i]
	
	start_com = np.array([[start_comx],[start_comy],[start_comz]])
	end_com = np.array([[end_comx,end_comy,end_comz]])
	
	filetypes.write_xyz('out',start_atoms)
	return startx,starty,startz,endx,endy,endz,start_com,end_com

def unweighted_COM(start_atoms,end_atoms):
	start_comx = sum(a.x for a in start_atoms)/len(start_atoms)
	start_comy = sum(a.y for a in start_atoms)/len(start_atoms)
	start_comz = sum(a.z for a in start_atoms)/len(start_atoms)
	
	end_comx = sum(a.x for a in end_atoms)/len(end_atoms)
	end_comy = sum(a.y for a in end_atoms)/len(end_atoms)
	end_comz = sum(a.z for a in end_atoms)/len(end_atoms)
	
	startx = [a.x-start_comx for a in start_atoms]
	starty = [a.y-start_comy for a in start_atoms]
	startz = [a.z-start_comz for a in start_atoms]

	endx = [a.x-end_comx for a in end_atoms]
	endy = [a.y-end_comy for a in end_atoms]
	endz = [a.z-end_comz for a in end_atoms]
	
	for i in range(len(start_atoms)):
		start_atoms[i].x = startx[i]
		start_atoms[i].y = starty[i]
		start_atoms[i].z = startz[i]

	start_com = np.array([[start_comx],[start_comy],[start_comz]])
	end_com = np.array([[end_comx,end_comy,end_comz]])
	filetypes.write_xyz('out',start_atoms)
	return startx,starty,startz,endx,endy,endz,start_com,end_com
	
def SVD(dot_size,Niter):
	lattice_start = []
	lattice_end = []
	dot1 = [1, 2, 4, 5, 7, 9, 10]
	dot4 = [4, 5, 6, 8]
	for run in dot4:
		# W is the optimal rotation matrix that minimizes the differences between the start and end coordinates
		# read xyz file and compare the start and end structures
		f = open('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size)+'.xyz')
		Natoms = int(f.readline())
		f.close()
		atoms = filetypes.parse_xyz('/fs/home/ja539/Desktop/C/jovana_test/HSEHsol_test_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size)+'.xyz')
		start_atoms = atoms[-Natoms:]
		energy, end_atoms = gaussian.parse_atoms('HSEHsol_testopt4_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size),check_convergence=False)
		#energy, end_atoms = gaussian.parse_atoms('/fs/home/jms875/Documents/nanocrystals/pb_oleate_hydrate/gaussian/dot2_m06_6.log',check_convergence=False)
		#energy, end_atoms = gaussian.parse_atoms('gaussian/HSEHsol_core1_dot6_2.log')
	
		startx,starty,startz,endx,endy,endz,start_com,end_com = unweighted_COM(start_atoms,end_atoms)
		A = np.array([startx,starty,startz])
		#print A.shape
		B = np.array([endx,endy,endz])
		BT = B.transpose()
		C = np.dot(A,BT)
		#print C.shape
		
		U, s, V = np.linalg.svd(C, full_matrices=True)
		#print U.shape, V.shape
		if np.linalg.det(V) < 0:
			for i in range(3):
				V[2][i] =  V[2][i]*-1.0
		#print U.shape, V.shape
		W = np.dot(V.transpose(),U.transpose())	
		A_new = np.dot(W,A)
		
		# piece together an xyz file that has the end and shifted start one after the other
		for i in range(len(start_atoms)):
			start_atoms[i].x = A_new[0][i]
			start_atoms[i].y = A_new[1][i]
			start_atoms[i].z = A_new[2][i]
			end_atoms[i].x = B[0][i]
			end_atoms[i].y = B[1][i]
			end_atoms[i].z = B[2][i]
		
		make_xyz_file('SVD_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size), end_atoms, 'w')
		make_xyz_file('SVD_iter'+str(Niter)+'_'+str(run)+'_'+str(dot_size), start_atoms, 'a')
		
		# calculate RMS translation
		D = 0.0
		for a,b in zip(start_atoms,end_atoms):
			D += ((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z))**(0.5)
		#procrustes_distance = (D/len(start_atoms))**0.5
		procrustes_distance = D/len(start_atoms)
		#print procrustes_distance
		
		# relative separations
		# for the core:
		coreD = 0.0
		for a1,b1 in zip(start_atoms[0:dot_size*2],end_atoms[0:dot_size*2]):
			for a2,b2 in zip(start_atoms[0:dot_size*2],end_atoms[0:dot_size*2]):
				coreD += abs(((a1.x-a2.x)*(a1.x-a2.x) + (a1.y-a2.y)*(a1.y-a2.y) + (a1.z-a2.z)*(a1.z-a2.z))**(0.5) - ((b1.x-b2.x)*(b1.x-b2.x) + (b1.y-b2.y)*(b1.y-b2.y) + (b1.z-b2.z)*(b1.z-b2.z))**(0.5))
		
		coreD = coreD/(dot_size*2)/(dot_size*2-1)
		print 'core:', coreD
		
		# for the ligands (average over all the ligands)
		ligD = [0.0]*((Natoms-2*dot_size)/10)
		for i in range((Natoms-2*dot_size)/10):
			for a1,b1 in zip(start_atoms[dot_size*2+10*i:dot_size*2+10*(i+1)],end_atoms[dot_size*2+10*i:dot_size*2+10*(i+1)]):
				for a2,b2 in zip(start_atoms[dot_size*2+10*i:dot_size*2+10*(i+1)],end_atoms[dot_size*2+10*i:dot_size*2+10*(i+1)]):
					ligD[i] += abs(((a1.x-a2.x)*(a1.x-a2.x) + (a1.y-a2.y)*(a1.y-a2.y) + (a1.z-a2.z)*(a1.z-a2.z))**(0.5) - ((b1.x-b2.x)*(b1.x-b2.x) + (b1.y-b2.y)*(b1.y-b2.y) + (b1.z-b2.z)*(b1.z-b2.z))**(0.5))
		
			ligD[i] = ligD[i]/10/9
		#print ligD
		ligDavg = (sum(a for a in ligD))/len(ligD)
		print 'lig:', ligDavg
		allD = 0.0
		for a1,b1 in zip(start_atoms,end_atoms):
			for a2,b2 in zip(start_atoms,end_atoms):
				allD += abs(((a1.x-a2.x)*(a1.x-a2.x) + (a1.y-a2.y)*(a1.y-a2.y) + (a1.z-a2.z)*(a1.z-a2.z))**(0.5) - ((b1.x-b2.x)*(b1.x-b2.x) + (b1.y-b2.y)*(b1.y-b2.y) + (b1.z-b2.z)*(b1.z-b2.z))**(0.5))
		
		allD = allD/Natoms/(Natoms-1)
		print 'all:', allD
		
		
		for a1,b1 in zip(start_atoms[0:dot_size],end_atoms[0:dot_size]):
			for a2,b2 in zip(start_atoms[dot_size:dot_size*2],end_atoms[dot_size:dot_size*2]):
				lattice_start += [((a1.x-a2.x)*(a1.x-a2.x) + (a1.y-a2.y)*(a1.y-a2.y) + (a1.z-a2.z)*(a1.z-a2.z))**(0.5)]
				lattice_end += [((b1.x-b2.x)*(b1.x-b2.x) + (b1.y-b2.y)*(b1.y-b2.y) + (b1.z-b2.z)*(b1.z-b2.z))**(0.5)]
		
	#print lattice_start
	#print lattice_end
	
SVD(4,4)
