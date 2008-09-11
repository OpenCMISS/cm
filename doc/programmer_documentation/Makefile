# define a makefile variable for the java 
#
JAVA = java
JFLAG = -jar
  
# Removes all .html files, so that the next make rebuilds them
#
clean: 
	$(RM) *.html

single:
	$(JAVA) $(JFLAG) lib/xalan.jar -in xml/toplevel.xml -out html/single.xhtml  -xsl nwalsh/xhtml/docbook.xsl

chunk:
	$(JAVA) $(JFLAG) lib/xalan.jar -in xml/toplevel.xml -out html/index.xhtml -xsl nwalsh/xhtml/chunk.xsl -param chunk.section.depth 0 -param section.autolabel 1
	
