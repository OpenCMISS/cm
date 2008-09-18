# define a makefile variable for the java 
#
JAVA = java
JFLAG = -jar
  
# Removes all .html files, so that the next make rebuilds them
#
clean: 
	$(RM) html/*.html

single:
ifdef IMPORT_DOXYGEN
	$(JAVA) $(JFLAG) lib/doxyToDoc.jar xml/library_commands.xml_template ../doxygen/xml
endif
	$(JAVA) $(JFLAG) lib/xalan.jar -in xml/toplevel.xml -out html/single.xhtml  -xsl customized_style/xhtml/customized-single.xsl

chunk:
ifdef IMPORT_DOXYGEN
	$(JAVA) $(JFLAG) lib/doxyToDoc.jar xml/library_commands.xml_template ../doxygen/xml
endif
	$(JAVA) $(JFLAG) lib/xalan.jar -in xml/toplevel.xml -out html/index.xhtml -xsl customized_style/xhtml/customized-chunk.xsl -param chunk.section.depth 0 -param section.autolabel 1
	
fo:
	$(JAVA) $(JFLAG) lib/xalan.jar -in xml/toplevel.xml -out pdf/programmer_doc.fo -xsl customized_style/fo/customized-fo.xsl -param section.autolabel 1 -param paper.type A4
	$(JAVA) $(JFLAG) lib/fop.jar pdf/programmer_doc.fo pdf/programmer_doc.pdf   
