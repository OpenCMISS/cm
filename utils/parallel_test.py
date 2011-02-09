#!/usr/bin/env python

"""
Run an OpenCMISS example in parallel with a profiling tool, then
parse the output
"""

import subprocess
import os,sys,errno
import shutil

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
        flatprocs=flatten(self.num_procs)
        if tool == "callgrind":
            results['totals'] = {}
            for func in function_list:
                results[func] = [0.0]*len(flatprocs)
            for (i,procs) in enumerate(flatprocs):
                test_outputdir=self.outputdir+os.sep+str(procs)
                callgrindfiles=[f for f in os.listdir(test_outputdir) if f.startswith('callgrind.out.')]
                for f in callgrindfiles:
                    #parse callgrind output with callgrind_annotate
                    subprocess.call('callgrind_annotate --inclusive=yes '+test_outputdir+os.sep+f+' > '+test_outputdir+os.sep+'annotate.'+f,shell=True)
                annotatedfiles=['annotate.'+f for f in callgrindfiles]
                #get average total number of instructions
                results['totals'][i]=self.__getcgtotal(test_outputdir,callgrindfiles)
                for func in function_list:
                    #get average for functions
                    results[func][i]=self.__getcgfunction(test_outputdir,annotatedfiles,callgrindfiles,func)
        elif tool == "massif":
            results['totals'] = {}
            for func in function_list:
                results[func] = [0.0]*len(flatprocs)
            for (i,procs) in enumerate(flatprocs):
                test_outputdir=self.outputdir+os.sep+str(procs)
                massiffiles=[f for f in os.listdir(test_outputdir) if f.startswith('massif.out.')]
                results['totals'][i]=self.__getmassiftotal(test_outputdir,massiffiles)
                for func in function_list:
                    results[func][i]=self.__getmassiffunction(test_outputdir,massiffiles,func)
        elif tool == "time":
            raise RuntimeError, "Not implemented"
        else:
            raise RuntimeError, "Invalid tool specified"

    def __getmassiftotal(self,dir,files):
        """Get the total peak memory usage from massif output"""
        values=[]
        for f in files:
            output=open(dir+os.sep+f,'r')
            for l in output.readlines():
                if l.find("mem_heap_B=")==0:
                    mem_usage=int(l.strip().split('=')[1])
                if l.find("heap_tree=peak")==0:
                    values.append(mem_usage)
                    break
        if len(values)>0:
            return float(sum(values))/float(len(values))
        else: return None

    def __getmassiffunction(self,dir,files,func):
        """Find the maximum memory used by a function from the massif output"""
        values=[]
        for f in files:
            output=open(dir+os.sep+f,'r')
            func_values=[]
            for l in output.readlines():
                if l.find(func) > -1:
                    func_values.append(int(l.split()[1]))
            if len(func_values)>0:
                values.append(max(func_values))
        if len(values)>0:
            return float(sum(values))/float(len(values))
        else: return None

    def __getcgtotal(self,dir,callgrindfiles):
        """Get the average total number of CPU instructions from callgrind output"""
        values=[]
        for f in callgrindfiles:
            cgoutput=open(dir+os.sep+f,'r')
            for l in cgoutput.readlines():
                if l.startswith('summary:'):
                    values.append(int(l.split()[1]))
                    break
            cgoutput.close()
        if len(values)>0:
            return float(sum(values))/float(len(values))
        else: return None

    def __getcgfunction(self,dir,annotatedfiles,callgrindfiles,func):
        """Get the average number of CPU instructions from callgrind output for a function"""
        values=[]
        for (a,c) in zip(annotatedfiles,callgrindfiles):
            value=0
            aoutput=open(dir+os.sep+a,'r')
            cgoutput=open(dir+os.sep+c,'r')
            # look in callgrind annotate output first
            for l in aoutput.readlines():
                if l.find(func) > -1:
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
                        if len(line) == 2 and line[0] in ('*','0'):
                            value = int(line[1])
                        break
                        # if it wasn't found then we'd have to sum up all the subroutines called from
                        # this subroutine, but we won't bother
                    else:
                        if l.find(func) > -1:
                            # skip 1 line to get second line after function name
                            step = 1
            if value:
                values.append(value)
            aoutput.close()
            cgoutput.close()
        if len(values)>0:
            return float(sum(values))/float(len(values))
        else: return None

    def write_results(self,outputfile):
        """Output csv file with results"""
        output=open(outputfile,'w')
        output.write(',')
        flatprocs=flatten(self.num_procs)
        for procs in flatprocs:
            output.write(str(procs)+',')
        output.write('\n')
        for tool in self.results.keys():
            output.write('"'+tool+'"\n')
            #output totals in first row
            if self.results[tool].has_key('totals'):
                func='totals'
                output.write('"'+func+'",')
                for i in range(len(flatprocs)):
                    output.write(str(self.results[tool][func][i])+',')
                output.write('\n')
            for func in self.results[tool].keys():
                #now output other functions
                if func != 'totals':
                    output.write('"'+func+'",')
                    for i in range(len(flatprocs)):
                        output.write(str(self.results[tool][func][i])+',')
                    output.write('\n')
            output.write('\n')
        output.close()

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
    profiler.runtest('massif')

    cg_funclist=["generated_mesh_create_finish",
        "decomposition_create_finish",
        "generated_mesh_geometric_parameters_calculate",
        "boundary_conditions_create_finish",
        "equations_set_assemble_static_linear_fem",
        "solver_linear_iterative_solve",
        "problem_solver_equations_solve",
        "cmiss_finalise",
        "problem_solver_equations_create_finish",
        "field_create_finish"]
    massif_funclist=[ "__field_routines_MOD_field_mappings_calculate", "__mesh_routines_MOD_mesh_topology_elements_adjacent_elements_calculate",
        "__mesh_routines_MOD_mesh_topology_nodes_calculate", "__mesh_routines_MOD_decomposition_topology_lines_calculate",
        "__mesh_routines_MOD_mesh_topology_elements_create_start", "__mesh_routines_MOD_decomp_topology_elem_adjacent_elem_calculate",
        "__mesh_routines_MOD_domain_topology_initialise_from_mesh", "__lists_MOD_list_initialise", "__mesh_routines_MOD_domain_mappings_nodes_dofs_calculate",
        "__solver_mapping_routines_MOD_solver_mapping_calculate" ]
    profiler.parse_results('callgrind',cg_funclist)
    profiler.parse_results('massif',massif_funclist)

    profiler.write_results(outputdir+os.sep+'output.csv')
