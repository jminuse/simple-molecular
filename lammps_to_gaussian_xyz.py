def convert_xyz(filename):
	f = open('lammps/'+filename+'.xyz')
	g = open(filename+'.xyz','w')
	lines = f.readlines()
	f.close()
	for line in lines:
		parts = line.split()
		if len(parts) == 4:
			if parts[0] == '3':
				parts[0] = 'S'
			elif parts[0] == '1' or parts[0] == '2':
				parts[0] = 'Pb'
			elif parts[0] == '5' or parts[0] == '4':
				parts[0] = 'O'
			elif parts[0] == '6' or parts[0] == '7' or parts[0] == '8'  or parts[0] == '9':
				parts[0] = 'C'
			else:
				parts[0] = 'H'
		
		for part in parts:
			g.write(part+' ')
		
		g.write('\n')
		
	g.close()

def convert_trunc_xyz(filename):
	f = open('lammps/'+filename+'.xyz')
	g = open(filename+'.xyz','w')
	lines = f.readlines()
	f.close()
	for line in lines:
		parts = line.split()
		if len(parts) == 4:
			if parts[0] == '3':
				parts[0] = 'S'
			elif parts[0] == '1' or parts[0] == '2':
				parts[0] = 'Pb'
			elif parts[0] == '5' or parts[0] == '4':
				parts[0] = 'O'
			elif parts[0] == '6' or parts[0] == '7':
				parts[0] = 'C'
			else:
				parts[0] = 'H'
		
		for part in parts:
			g.write(part+' ')
		
		g.write('\n')
		
	g.close()

def combine_xyz(filename_prefix):
	good_dot2 = [[1,2,3,4,5,6,8],[1,2,6,8,9],[1,3,4,5,7,8,9],[1,2,3,4,6,7,8],[1,3,4,5,6,8],[1,2,4,5,6,7,8,9],[1,2,3,4,5,6,7,8],[1,2,3,4,6,7,8,9],[1,2,4,5,6,7,8,9],[1,2,3,6,7,8,9]]
	dot1 = 1
	for seed in range(1,10):
		for dot2 in good_dot2[seed-1]:
			g = open(filename_prefix+'s'+str(seed)+'_spc_'+str(dot1)+'_'+str(dot2)+'.xyz','w')
			for step in range(1,5):
				f = open('lammps/'+filename_prefix+'s'+str(seed)+'_spc'+str(step)+'_'+str(dot1)+'_'+str(dot2)+'.xyz')
				lines = f.readlines()
				f.close()
				for line in lines:
					parts = line.split()
					if len(parts) == 4:
						if parts[0] == '3':
							parts[0] = 'S'
						elif parts[0] == '1' or parts[0] == '2':
							parts[0] = 'Pb'
						elif parts[0] == '5' or parts[0] == '4':
							parts[0] = 'O'
						elif parts[0] == '6' or parts[0] == '7' or parts[0] == '8'  or parts[0] == '9':
							parts[0] = 'C'
						else:
							parts[0] = 'H'
		
					for part in parts:
						g.write(part+' ')
		
					g.write('\n')
		
			g.close()


#convert_xyz('HSEH_dot3_new_20')
#convert_xyz('HSEHsol2_dot10')
def run_convert_gyr():
	for i in range(11,16):
		convert_xyz('gyr_s'+str(i)+'_1_1')
	
def run_convert_gyr_re():
	for i in range(1,5):
		convert_xyz('gyr_re_s'+str(i)+'_1_1')

def run_convert_runabf():
	for i in range(1,5):
		convert_xyz('gyr_runabf_s'+str(i)+'_1_1')
	
def run_convert_runabf2():
	for i in range(1,5):
		convert_xyz('gyr_runabf2_s'+str(i)+'_1_1')
	
def run_convert_p_abf():
	for dot2 in range(2,10):
		for seed in range(1,5):
			for part in range(1,5):
				convert_xyz('gyr_p'+str(part)+'_s'+str(seed)+'_1_'+str(dot2))

def run_convert_biasabf():
	for seed in range(1,5):
		convert_xyz('gyr_biasabf_s%d_1_1' % seed)
		
def run_convert_joinabf():
	for seed in range(1,5):
		convert_xyz('gyr_joinabf_s%d_1_1' % seed)

def run_convert_joinabf2():
	for seed in range(1,5):
		convert_xyz('gyr_joinabf2_s%d_1_1' % seed)

def run_convert_joinabf2_25():
	for seed in range(1,5):
		convert_xyz('gyr_joinabf2_25_s%d_1_1' % seed)

def run_convert_joinabf2_5():
	for seed in range(1,5):
		convert_xyz('gyr_joinabf2_5_s%d_1_1' % seed)

def run_convert_joinabf3():
	for seed in range(1,5):
		convert_xyz('gyr_joinabf3_s%d_1_1' % seed)

def run_convert_meta_trunc():
	for i in range(4,5):
		convert_trunc_xyz('gyr_trunc_s'+str(i)+'_1_1')

def run_convert_meta_trunc_abf():
	for i in range(1,11):
		convert_trunc_xyz('gyr_trunc_abf_s'+str(i)+'_1_1')

