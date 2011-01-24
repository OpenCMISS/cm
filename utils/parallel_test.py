#!/usr/bin/env python

"""
Run an OpenCMISS example in parallel with a profiling tool, then
parse the output
"""

import subprocess
import os,sys,errno
import shutil
import re

def mkdire(path):
    """Make directory, including parents, and ignore error if it already exists"""
    try:
        os.makedirs(path)
    except OSError, err:
        if err.errno == errno.EEXIST:
            pass
        else: raise

def flatten(a):
    """Flatten a list of lists"""
    return [item for sublist in a for item in sublist]

class CmProfile():
    def __init__(self,num_procs,outputdir,source,executable,executable_options=""):
        self.num_procs=num_procs # must be a list of lists
        #eg num_procs=[[1,2,4,8],[16]] means run jobs with 1,2,4,8 processors in parallel concurrently,
        #then run with 16 processros once those jobs are finished
        self.outputdir=outputdir
        self.source=source
        self.executable=executable
        self.executable_options=executable_options

        self.massifopts='--threshold=0.1 --max-snapshots=200'
        self.callgrindopts=''
        self.timeopts=''

        self.results={}

    def savesettings(self):
        #copy the source file and save settings to keep a record of the test configuration
        settingsdir=self.outputdir+os.sep+'test_settings'
        mkdire(settingsdir)
        shutil.copy(self.source,settingsdir+'/source.f90')
        settingsfile=open(settingsdir+os.sep+'test_settings.txt','w')
        settingsfile.write('Source file: '+self.source+'\n')
        settingsfile.write('Executable: '+self.executable+'\n')
        settingsfile.write('Executable options: '+self.executable_options+'\n')
        settingsfile.write('Processors: '+str(self.num_procs)+'\n')
        settingsfile.close()

    def runtest(self,tool):
        if not os.path.isfile(self.executable):
            raise RuntimeError, "Executable file '%s' not found" % self.executable
        #check mpd is running and run if not
        ret = subprocess.call('mpdringtest',shell=True)
        if ret != 0:
            raise RuntimeError, "Error running mpd"
        for procs_list in self.num_procs:
            processes=[]
            for procs in procs_list:
                test_outputdir=self.outputdir+os.sep+str(procs)
                mkdire(test_outputdir)
                if tool == "callgrind":
                    cmd = 'mpiexec -n %d valgrind --tool=callgrind --log-file="%s/valgrind.out.%%p" ' \
                        '--callgrind-out-file="%s/callgrind.out.%%p" %s %s %s' %  \
                        (procs,test_outputdir,test_outputdir,self.callgrindopts,self.executable,self.executable_options)
                    processes.append(subprocess.Popen(cmd,shell=True))
                elif tool == "massif":
                    cmd = 'mpiexec -n %d valgrind --tool=massif --log-file="%s/valgrind.out.%%p" ' \
                        '--massif-out-file="%s/massif.out.%%p" %s %s %s' % \
                        (procs,test_outputdir,test_outputdir,self.massifopts,self.executable,self.executable_options)
                    processes.append(subprocess.Popen(cmd,shell=True))
                elif tool == "time":
                    # %e %U %S %M %x = real time, user time, system/kernel time, maximum resident memory, return value
                    # all times in seconds, memory in kb
                    cmd = 'mpiexec -n %d  /usr/bin/time -a -o %s/time_%d.out -f "%%e %%U %%S %%M %%x" %s %s %s' % \
                        (procs,test_outputdir,procs,self.timeopts,self.executable,self.executable_options)
                    processes.append(subprocess.Popen(cmd,shell=True))
                else:
                    raise RuntimeError, "Invalid tool specified"
            for process in processes:
                #wait for all jobs to finish before starting the next group of jobs
                process.wait()

    def parse_results(self,tool,function_list=[]):
        self.results[tool]={}
        results=self.results[tool]
        if tool == "callgrind":
            results['totals'] = {}
            for func in function_list:
                results[func] = {}
            for procs in flatten(self.num_procs):
                test_outputdir=self.outputdir+os.sep+str(procs)
                callgrindfiles=[f for f in os.listdir(test_outputdir) if f.startswith('callgrind.out.')]
                for f in callgrindfiles:
                    #parse callgrind output with callgrind_annotate
                    subprocess.call('callgrind_annotate --inclusive=yes '+test_outputdir+os.sep+f+' > '+test_outputdir+os.sep+'annotate.'+f,shell=True)
                annotatedfiles=['annotate.'+f for f in callgrindfiles]
                #get average total number of instructions
                results['totals'][procs]=self.__getcgtotal(test_outputdir,callgrindfiles)
                for func in function_list:
                    #get average for functions
                    results[func][procs]=self.__getcgfunction(test_outputdir,annotatedfiles,callgrindfiles,func)
        elif tool == "massif":
            raise RuntimeError, "Not implemented"
        elif tool == "time":
            raise RuntimeError, "Not implemented"
        else:
            raise RuntimeError, "Invalid tool specified"

    def __getcgtotal(self,dir,callgrindfiles):
        values=[]
        for f in callgrindfiles:
            cgoutput=open(dir+os.sep+f,'r')
            for l in cgoutput.readlines():
                if l.startswith('summary:'):
                    values.append(int(l.split()[1]))
                    break
            cgoutput.close()
        if len(values)>0:
            return sum(values)/len(values)
        else: return None

    def __getcgfunction(self,dir,annotatedfiles,callgrindfiles,func):
        values=[]
        for (a,c) in zip(annotatedfiles,callgrindfiles):
            value=0
            aoutput=open(dir+os.sep+a,'r')
            cgoutput=open(dir+os.sep+c,'r')
            # look in callgrind annotate output first
            for l in aoutput.readlines():
                m = re.search(func,l)
                if m:
                    value = int(l.split()[0].replace(',',''))
                    break
            if not value:
                # read from callgrind output
                step = -1
                for l in cgoutput.readlines():
                    if step>0:
                        step-=1
                    elif step==0:
                        line = l.split()
                        if len(line) == 2:
                            value = int(line[1])
                        break
                        # if it wasn't found then we'd have to sum up all the subroutines called from
                        # this subroutine, but we won't bother
                    else:
                        m = re.search(func,l)
                        if m:
                            # skip 1 line to get second line after function name
                            step = 1
            if value:
                values.append(value)
            aoutput.close()
            cgoutput.close()
        if len(values)>0:
            return sum(values)/len(values)
        else: return None

    def write_results(self,outputfile):
        """Output csv file with results"""
        output=open(outputfile,'w')
        output.write(',')
        for procs in flatten(self.num_procs):
            output.write(str(procs)+',')
        output.write('\n')
        for tool in self.results.keys():
            output.write('"'+tool+'"\n')
            for func in self.results[tool].keys():
                output.write('"'+func+'",')
                for procs in flatten(self.num_procs):
                    output.write(str(self.results[tool][func][procs])+',')
                output.write('\n')
            output.write('\n')
        output.close()

    def plot_results(self,plotdir):
        """Plot results and save plots as pngs in plotdir"""
        import matplotlib as plot
        raise RuntimeError, "Not implemented"

if __name__ == "__main__":
    #example usage
    num_procs=[[1,2,],[4]]
    options='10 10 10 1'
    outputdir='laplace_test'
    exampledir=os.environ['OPENCMISS_ROOT']+'/cm/examples/ClassicalField/Laplace/ParallelLaplace/'
    source=exampledir+'src/ParallelLaplaceExample.f90'
    executable=exampledir+'bin/x86_64-linux/mpich2/gnu_4.4/ParallelLaplaceExample'

    profiler = CmProfile(num_procs,outputdir,source,executable,options)
    profiler.savesettings()
    profiler.runtest('callgrind')

    funclist=["generated_mesh_create_finish",
        "decomposition_create_finish",
        "generated_mesh_geometric_parameters_calculate",
        "boundary_conditions_create_finish",
        "equations_set_assemble_static_linear_fem",
        "solver_linear_iterative_solve",
        "problem_solver_equations_solve",
        "cmiss_finalise",
        "problem_solver_equations_create_finish",
        "field_create_finish"]
    profiler.parse_results('callgrind',funclist)
    profiler.write_results(outputdir+os.sep+'output.csv')
