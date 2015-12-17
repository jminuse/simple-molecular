def convert_log(param_file):
# param file denotes where to read off the bounds for the params
	f = open(param_file+'.txt')
	g = open(param_file+'_converted.txt','w')
	lines = f.readlines()
	for i in range(len(params)):
		line = lines[i]
		line = line.split()
		g.write(str(params[i])+' '+line[1]+' '+line[2]+'\n')
		
	g.close()

	convert_log('HSEH_test_out')
	
def convert_bounds(param_file):
	# outputs bounds in a column of text
	f = open(param_file+'.txt')
	g = open(param_file+'_converted.txt','w')
	lines = f.readlines()
	for i in range(len(lines)):
		line = lines[i]
		line = line.split()
		g.write(line[1]+'\n')
	g.write('\n')
	for i in range(len(lines)):
		line = lines[i]
		line = line.split()
		g.write(line[2]+'\n')
		
	g.close()

	convert_bounds('HSEH_test3_1_out2')
