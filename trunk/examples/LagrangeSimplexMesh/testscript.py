import hashlib
import sys
import commands
import string

def md5sum(filename):
  myfile = open(filename)
  m = hashlib.md5()
  
  for aline in myfile :
    m.update(string.join(string.split(aline)))
  myfile.close()
  return m.hexdigest()

def checkFile(filename, compiler):
  expectedfile = "expected_files/"+filename+"_"+compiler
  md5actual =  md5sum(filename)
  md5expected =  md5sum(expectedfile)
  if md5actual == md5expected :
    print "LagrangeSimplexMesh - " + filename + " check succeeded."
    return 0
  else :
    print "LagrangeSimplexMesh - " + filename + " check FAILED!."
    return 1
  

if len(sys.argv)!=2 :
  print "Invalid number of arguments"
  sys.exit(1)
name = sys.argv[1]
print commands.getoutput(name)
namelist = string.split(name,'-')
compiler = namelist[len(namelist)-1]

if checkFile( "LagrangeSimplexMeshExample.part0.exelem", compiler ) == 1:
  sys.exit(1)

if checkFile( "LagrangeSimplexMeshExample.part0.exnode", compiler ) == 1:
  sys.exit(1)
