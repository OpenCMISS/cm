import FunctionRead

FunctionRead.main()

f_functions = open('opencmiss_f_functions', 'r')
c_functions = open('opencmiss_c_functions', 'r')
c_function = c_functions.readlines()
c_functions.close()

h_functions = open('opencmiss_h_functions', 'r')
h_function = h_functions.readlines()
h_functions.close()
missingfunctions_c = open('missing_functions_c', 'w')
missingfunctions_h = open('missing_functions_h', 'w')

countc = 0
counth = 0

string_f = f_functions.readline()

while string_f != '':
	matchC = False
	for string_c in c_function :
		if (string_f==string_c)	:
         		matchC=True
	if (not matchC):
		print 'Missing function in opencmiss_c.f90: ', string_f
		missingfunctions_c.write(string_f)
		countc = countc + 1
	matchH = False
	for string_h in h_function :
		if (string_f==string_h) :	
         		matchH=True
	if (not matchH):
		print 'Missing function in opencmiss.h: ', string_f		
		missingfunctions_h.write(string_f)
		counth = counth + 1

	string_f = f_functions.readline()

if countc ==0:
	print 'No functions missing in opencmiss_c.f90'
else:
	print countc, ' functions missing from opencmiss_c.f90'
if counth == 0:
	print 'No functions missing in opencmiss.h' 
else:
	print counth, ' functions missing from opencmiss.h'

missingfunctions_c.close()
missingfunctions_h.close()

f_functions.close()

