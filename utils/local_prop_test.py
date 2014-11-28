#!/usr/bin/env python
import sys,os
sys.path.append(os.environ['OPENCMISS_ROOT']+'/cm/examples')
from noseMain import *

"""
Use the noseMain.py testing file to run the tests in a .prop file locally
Doesn't build the example or check building
"""

def test_prop_file(prop_file_path):
  global compiler
  root = os.getcwd()
  system = os.uname()[0].lower()
  arch = os.uname()[4]
  propFile= file(prop_file_path, "r")
  properties = dict()
  load_prop(propFile,properties)
  if '42TestingPointsPATH' in properties :
    testingPointsPath = os.environ['OPENCMISS_ROOT']+"/cm/examples/"+properties['42TestingPointsPATH']
  else:
    testingPointsPath = ""
  testpoints = properties['TestingPoint']
  for testpoint in testpoints :
    if (testpoint[0].find("${")!=-1) :
      testingPointPath = properties[testpoint[0][2:testpoint[0].find("}")]]+testpoint[0][testpoint[0].find("}")+1:]
    else:
      testingPointPath = testpoint[0]
    os.chdir(testingPointsPath + testingPointPath)
    if len(testpoint)<=2 :            
      yield check_run, 'run', os.getcwd(), system, arch, compiler, root, testpoint[1], testingPointPath, False
    else :
      yield check_run, 'run', os.getcwd(), system, arch, compiler, root, testpoint[1], testingPointPath
  
if __name__ == '__main__':
  if len(sys.argv) != 2:
    sys.stderr.write('Usage '+sys.argv[0]+' prop_file\n')
    exit(1)
  for check in test_prop_file(sys.argv[1]):
      print ' '.join([str(c) for c in check[1:]])
      check[0](*check[1:])
