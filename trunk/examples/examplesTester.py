# -*- coding: utf-8 -*-
import os, sys, subprocess

cwd = os.getcwd();
logDir = cwd + "/../../build/logs";
rootUrl = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/"

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];
os.putenv('HOME', '/home/autotest')
mpidir = cwd+'/../../../opencmissextras/cm/external/x86_64-linux-debug/mpich2/'+compiler+'/bin'
os.system('python ' + mpidir + '/mpd.py &')
f = open(logDir+'/successBuilds',"r")
successbuilds = f.read()
f.close()
os.remove(logDir+'/successBuilds')

def testExample(id, path, nodes, input=None, args=None) :
   global compiler,logDir,successbuilds;
   if (successbuilds.find(path)!=-1) :
     newDir = logDir
     for folder in path.split('/') :
       newDir = newDir + '/' + folder
       if not os.path.isdir(newDir):
         os.mkdir(newDir)
     os.chdir(path)
     if os.path.exists(newDir + "/test"+id+"-" + compiler) :
       os.remove(newDir + "/test"+id+"-" + compiler)
     execPath='bin/x86_64-linux/mpich2/'+compiler+'/'+path.rpartition('/')[2]+'Example-debug'
     if nodes == '1' :
       if input != None :
         inputPipe = subprocess.Popen(["echo", input], stdout=subprocess.PIPE)
         f = open(newDir + "/test" + id + "-" + compiler,"w")
         execCommand = subprocess.Popen([execPath], stdin=inputPipe.stdout, stdout=f,stderr=f)
         f.close()
         err = os.waitpid(execCommand.pid, 0)[1]
       elif args==None :
         err=os.system(execPath +" > " + newDir + "/test" + id + "-" + compiler + " 2>&1")
       else :
         err=os.system(execPath + ' ' + args +" > " + newDir + "/test" + id + "-" + compiler + " 2>&1")
     else :
       if input != None :
         inputPipe = subprocess.Popen(["echo", input], stdout=subprocess.PIPE)
         f = open(newDir + "/test" + id + "-" + compiler,"w")
         execCommand = subprocess.Popen(["mpiexec","-n",nodes,execPath], stdin=inputPipe.stdout, stdout=f,stderr=subprocess.PIPE)
         f.close()
         err = os.waitpid(execCommand.pid, 0)[1]
       elif args==None :
         err=os.system('python ' + mpidir + '/mpiexec.py -n ' + nodes + ' ' + execPath +" > " + newDir + "/test" + id + "-" + compiler + " 2>&1")
       else :
         err=os.system('python ' + mpidir + '/mpiexec.py -n ' + nodes + " " + execPath + ' ' + args+" > " + newDir + "/test" + id + "-" + compiler + " 2>&1")
     if not os.path.exists(execPath) :
       err=-1
     if err==0 :
       print "Testing %s%s: <a class='success' href='%slogs_x86_64-linux/%s/test%s-%s'>success</a><br>" %(path,id,rootUrl,path,id,compiler)
     else :
       print "Testing %s%s: <a class='fail' href='%slogs_x86_64-linux/%s/test%s-%s'>failed</a><br>" %(path,id,rootUrl,path,id,compiler)
     os.chdir(cwd)
   else :
     print "Testing %s%s: <a class='fail'>failed</a> due to build failure<br>" %(path,id)
   return;

testExample(id='1', path="ClassicalField/AnalyticLaplace", nodes='1')
testExample(id='1', path="ClassicalField/AdvectionDiffusion", nodes='1')
testExample(id='1', path="ClassicalField/Diffusion", nodes='1')
testExample(id='1', path="ClassicalField/DiffusionConstantSource", nodes='1')
#testExample(id='1', path="ClassicalField/Helmholtz", nodes='1',input='4\n4\n0\n1')
#testExample(id='1', path="ClassicalField/Laplace", nodes='1',args='4 4 0 1') 
#testExample(id='2', path="ClassicalField/Laplace", nodes='2',args='4 4 0 2')
testExample(id='1', path="ClassicalField/NonlinearPoisson",nodes='1',input='4\n4\n0\n1')
testExample(id='1', path="ClassicalField/NewLaplace", nodes='1',args='4 4 0 1') 
testExample(id='2', path="ClassicalField/NumberLaplace", nodes='1',args='4 4 0 1')

#testExample(id='1',path="Bioelectrics/Monodomain",nodes='1',input='4\n4\n0\n1')
  
testExample(id='1',path="FluidMechanics/Stokes/ALE",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/Stokes/Static",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/Stokes/Dynamic",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/NavierStokes/ALE",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/NavierStokes/Static",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/NavierStokes/Dynamic",nodes='1',input='\n')


#testExample(id='1',path="FluidMechanics/Darcy/ConvergenceStudy",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="FluidMechanics/Darcy/FiveSpotProblem",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="FluidMechanics/Darcy/VenousCompartment",nodes='1',input='4\n4\n0\n1')

#testExample(id='1',path="FiniteElasticity/UniAxialExtension",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="FiniteElasticity/TwoElementTriLinear",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="FiniteElasticity/MixedBoundaryConditions",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="FiniteElasticity/TriCubicAxialExtension",nodes='1',input='4\n4\n0\n1')

#testExample(id='1',path="LinearElasticity/2DAnalytic1",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/2DPlaneStressLagrangeBasis",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/3DAnalytic1",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/3DLagrangeBasis",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/3DCubicHermiteBasis",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/3DLagrangeBasisAnisotropicFibreField",nodes='1',input='4\n4\n0\n1')


# TODO Group them
#testExample(id='1',path="LagrangeSimplexMesh",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="cellml",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="define-geometry-and-export",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="MoreComplexMesh",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="simple-field-manipulation-direct-access",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="SimplexMesh",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="TwoRegions",nodes='1',input='4\n4\n0\n1')

os.system('python ' + mpidir + '/mpdallexit.py')

