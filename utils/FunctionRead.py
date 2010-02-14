opencmiss_f = open('opencmiss.f90', 'r')
opencmiss_c = open('opencmiss_c.f90', 'r')
opencmiss_h = open('opencmiss.h', 'r')
f_functions = open('opencmiss_f_functions', 'w+')
c_functions = open('opencmiss_c_functions', 'w+')
h_functions = open('opencmiss_h_functions', 'w+')


def Read(opencmiss, f, startstring_ind, endstring_ind, endfunc_ind, startleng, endfunc_num, endfunc_ptr):

	print startstring_ind

	string_c = opencmiss.readline()
	endoffile = string_c.find(endfunc_ind)

	while endoffile == -1:

		string_c = opencmiss.readline()

		startstring = string_c.find(startstring_ind)
		endstring = string_c.find(endstring_ind)

		readstring = []

		if endstring != -1:
			if startstring != -1:
				for i in range(startstring+startleng, endstring):
					readstring.extend(string_c[i])
			joinstring = "".join(readstring)
			joinstring = joinstring.rstrip()
			if joinstring.find('C') != -1:
				if joinstring.endswith(endfunc_num):
					joinstring = joinstring.rstrip(endfunc_num)
				if joinstring.endswith(endfunc_ptr):
					joinstring = joinstring.rstrip(endfunc_ptr)
				if joinstring.endswith('C'):
					joinstring = joinstring.rstrip('C')
				print joinstring
				f.write(joinstring + '\n')
		endoffile = string_c.find(endfunc_ind)
	f.write('endoffile')

Read(opencmiss_f, f_functions, 'SUBROUTINE CMISS', '(', 'END MODULE', 11, 'Number', 'Obj')
Read(opencmiss_c, c_functions, 'FUNCTION CMISS', '(', 'END MODULE', 9, 'Num', 'Ptr')
Read(opencmiss_h, h_functions, 'CMISSError CMISS', '(', 'EndOfFile', 11, 'Num', '')

opencmiss_c.close()
c_functions.close()
opencmiss_f.close()
f_functions.close()
opencmiss_h.close()
h_functions.close()
