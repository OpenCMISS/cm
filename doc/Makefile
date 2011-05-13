# Makefile to build documentations
#----------------------------------------------------------------------------------------------------------------------------------

clean:
	find ../../cm -type d -exec chmod 755 {} \;
	rm -rf doxy/latexhtml doxy/programmer doxy/user doxy/OpenCMISSNotes.pdf

relocate:
	mv notes/OpenCMISSNotes.pdf doxy
	mv notes/OpenCMISSNotes doxy/latexhtml
	mv doxygen/Programmer/html doxy/programmer
	mv doxygen/User/html doxy/user

install:
	cd doxygen; \
	doxygen Doxyfile; \
	doxygen Doxyfile_User
	cd notes; \
	./Latex_make.sh pdf; \
	./Latex_make.sh html


  
