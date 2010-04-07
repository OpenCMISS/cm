def suffixOut(str1, suffixes) :
	for suffix in suffixes :
        	if str1.endswith(suffix):
			return suffix
	return ''

def Read(opencmiss, f, startstring_ind, endstring_ind, startleng, endfunc):

	string_c = opencmiss.readline()
        string_list = []

	while string_c != '':

		string_c = opencmiss.readline()

		if not string_c.startswith('!') :
			startstring = string_c.find(startstring_ind)
			endstring = string_c.find(endstring_ind)

			readstring = []

			if endstring != -1:
				if startstring != -1:
					for i in range(startstring+startleng, endstring):
						readstring.extend(string_c[i])
				joinstring = "".join(readstring)
				joinstring = joinstring.rstrip()
				suffix = suffixOut(joinstring, endfunc)
				while suffix != '' :
					joinstring = joinstring[:-len(suffix)]
					suffix = suffixOut(joinstring, endfunc)
				if string_list.count(joinstring)==0 and joinstring!='' :
                        		string_list.append(joinstring)
					f.write(joinstring + '\n')

def main() :
	opencmiss_f = open('../src/opencmiss.f90', 'r')
	opencmiss_c = open('../src/opencmiss_c.f90', 'r')
	opencmiss_h = open('../src/opencmiss.h', 'r')
	f_functions = open('opencmiss_f_functions', 'w+')
	c_functions = open('opencmiss_c_functions', 'w+')
	h_functions = open('opencmiss_h_functions', 'w+')

	Read(opencmiss_f, f_functions, 'SUBROUTINE CMISS', '(', 11, ['Number','C', 'Obj','VS','0','1'])
	Read(opencmiss_c, c_functions, 'FUNCTION CMISS', '(', 9, ['Num', 'Ptr', 'C'])
	Read(opencmiss_h, h_functions, 'CMISSError CMISS', '(', 11, ['Num','C'])

	opencmiss_c.close()
	c_functions.close()
	opencmiss_f.close()
	f_functions.close()
	opencmiss_h.close()
	h_functions.close()
