#test

f_functions = open('opencmiss_f_functions', 'r')
c_functions = open('opencmiss_c_functions', 'r')
h_functions = open('opencmiss_h_functions', 'r')
missingfunctions_c = open('missing_functions_c', 'w')
missingfunctions_h = open('missing_functions_h', 'w')
functionchecks_c = open('functionchecks_c', 'w')

endoffile = -1
endoffilec = -1
endoffileh = -1

countc = 0
counth = 0
comparestringc = 1
comparestringh = 1

while endoffile == -1:

	string_f = f_functions.readline()
	endoffile = string_f.find('endoffile')
	string_c = c_functions.readline()
	endoffilec = string_c.find('endoffile')
	string_h = h_functions.readline()
	endoffileh = string_h.find('endoffile')

	while string_f != string_c and endoffilec == -1:
		endoffilec = string_c.find('endoffile')
		string_c = c_functions.readline()
	c_functions.seek(0)
	string_c = c_functions.readline()
	if endoffilec != -1:
		missingfunctions_c.write(string_f + '\n')
		countc = countc + 1

	while string_f != string_h and endoffileh == -1:
		endoffileh = string_h.find('endoffile')
		string_h = h_functions.readline()
	h_functions.seek(0)
	string_h = h_functions.readline()
	
	if string_f == string_h:
			print 'Function ' + string_f + 'matched to ' + string_h
	if endoffileh != -1:
		missingfunctions_h.write(string_f + '\n')
		counth = counth + 1

	endoffile = string_f.find('endoffile')
	string_f = f_functions.readline()

if countc ==0:
	print 'No functions missing in opencmiss_c.f90'
else:
	print countc, ' functions missing from opencmiss_c.f90'
if counth == 0:
	print 'No functions missing in opencmiss.h' 
else:
	print countc, ' functions missing from opencmiss.h'

missingfunctions_c.close()
missingfunctions_h.close()
c_functions.close()
f_functions.close()
h_functions.close()
