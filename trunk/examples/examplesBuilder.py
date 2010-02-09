# -*- coding: utf-8 -*-
import os, sys

cwd = os.getcwd();
logDir = cwd + "/../../build/logs";
rootUrl = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/"

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];
f = open(logDir+'/successBuilds',"w")

def buildExample(path) :
   global compiler,logDir,libSuccess;
   if (libSuccess!=-1) :
     newDir = logDir
     for folder in path.split('/') :
       newDir = newDir + '/' + folder
       if not os.path.isdir(newDir):
         os.mkdir(newDir)
     os.chdir(path)
     if os.path.exists(newDir + "/build-" + compiler) :
       os.remove(newDir + "/build-" + compiler)
     err=os.system("make COMPILER=" + compiler + " > " + newDir + "/build-" + compiler +" 2>&1")
     if err==0 :
       f.write(path+'\n')
       print "Building %s: <a class='success' href='%slogs_x86_64-linux/%s/build-%s'>success</a><br>" %(path,rootUrl,path,compiler)
     else :
       print "Building %s: <a class='fail' href='%slogs_x86_64-linux/%s/build-%s'>failed</a><br>" %(path,rootUrl,path,compiler)
     os.chdir(cwd)
   else :
     print "Building %s: <a class='fail'>failed</a> due to library build failure<br>" %(path)
   return;

def buildLibrary() :
   global compiler,logDir;
   newDir = logDir
   os.chdir('..')
   if os.path.exists(newDir + "/build-" + compiler) :
     os.remove(newDir + "/build-" + compiler)
   err=os.system("make COMPILER=" + compiler + " > " + newDir + "/build-" + compiler +" 2>&1")
   if err==0 :
     print "Building OpenCMISS Library: <a class='success' href='%slogs_x86_64-linux/build-%s'>success</a><br>" %(rootUrl,compiler)
     os.chdir(cwd)
     return 0;
   else :
     print "Building OpenCMISS Library: <a class='fail' href='%slogs_x86_64-linux/build-%s'>failed</a><br>" %(rootUrl,compiler)
     os.chdir(cwd)
     return -1;
   
   
libSuccess = buildLibrary()

buildExample("ClassicalField/AnalyticLaplace")
buildExample("ClassicalField/AdvectionDiffusion")
buildExample("ClassicalField/Diffusion")
buildExample("ClassicalField/DiffusionConstantSource")
#buildExample("ClassicalField/Helmholtz")
#buildExample("ClassicalField/Laplace")
buildExample("ClassicalField/NonlinearPoisson")
buildExample("ClassicalField/NewLaplace")
buildExample("ClassicalField/NumberLaplace")

#buildExample("Bioelectrics/Monodomain")

buildExample("FluidMechanics/Stokes/ALE")
buildExample("FluidMechanics/Stokes/Static")
buildExample("FluidMechanics/Stokes/Dynamic")

buildExample("FluidMechanics/NavierStokes/ALE")
buildExample("FluidMechanics/NavierStokes/Static")
buildExample("FluidMechanics/NavierStokes/Dynamic")


#buildExample("FluidMechanics/Darcy/ConvergenceStudy")
#buildExample("FluidMechanics/Darcy/FiveSpotProblem")
#buildExample("FluidMechanics/Darcy/VenousCompartment")

#buildExample("FiniteElasticity/UniAxialExtension")
#buildExample("FiniteElasticity/TwoElementTriLinear")
#buildExample("FiniteElasticity/MixedBoundaryConditions")
#buildExample("FiniteElasticity/TriCubicAxialExtension")

#buildExample("LinearElasticity/2DAnalytic1")
#buildExample("LinearElasticity/2DPlaneStressLagrangeBasis")
#buildExample("LinearElasticity/3DAnalytic1")
#buildExample("LinearElasticity/3DLagrangeBasis")
#buildExample("LinearElasticity/3DCubicHermiteBasis")
#buildExample("LinearElasticity/3DLagrangeBasisAnisotropicFibreField")


# TODO Group them
#buildExample("LagrangeSimplexMesh")
#buildExample("cellml")
#buildExample("define-geometry-and-export")
#buildExample("MoreComplexMesh")
#buildExample("simple-field-manipulation-direct-access")
#buildExample("SimplexMesh")
#buildExample("TwoRegions")

f.close()
