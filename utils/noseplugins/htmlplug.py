"""This is a very basic example of a plugin that controls all test
output. In this case, it formats the output as ugly unstyled html.

Upgrading this plugin into one that uses a template and css to produce
nice-looking, easily-modifiable html output is left as an exercise for
the reader who would like to see his or her name in the nose AUTHORS file.
"""
import traceback
from nose.plugins import Plugin
from datetime import date

import inspect,os

class ResultTree:
  isPassed = True
  parent = None
  childTrees = []

  def __init__(self,parent=None) :
    self.childTrees = []
    self.isPassed = True
    if parent!=None :
      self.parent=parent
      self.parent.childTrees.append(self)
  
  def isPass(self):
    if len(self.childTrees)==0 :
      return self.isPassed
    else :
      for childTree in self.childTrees:
        if not childTree.isPass() :
          return False
      return True
      
  def setPass(self, isPassed) :
    self.isPassed = isPassed


class HtmlOutput(Plugin):
    """Output test results as ugly, unstyled html.
    """
    
    name = 'html-output'
    current = ResultTree()
    testLevelsInner = []
    testLevelsOuter = []
    isStart = True
    buildbotUrl = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/"
    
    def __init__(self):
        super(HtmlOutput, self).__init__()
        self.header=['<html><head>',
                      '<title>Test output</title>',
                      '<SCRIPT LANGUAGE="JavaScript" SRC="http://autotest.bioeng.auckland.ac.nz/opencmiss-build/tree.js"> </SCRIPT>',
                      '<link href="http://autotest.bioeng.auckland.ac.nz/opencmiss-build/tree.css" rel="stylesheet" />',
                      '</head><body>',
                      '<h1>OpenCMISS Nighly Testing Results</h1>']
        self.html=['</div>','<ul id="tree1" class="mktree">']

    def insertLog(self, test):
        logPath=""
        historyPath=""
        if str(test).find('test_build_library')!=-1 :
          (compiler_version, arch, system) = list(test.test.arg)[:3]
          logPath = "%slogs_%s-%s/nose_library_build_%s%s" %(self.buildbotUrl, arch, system, compiler_version, str(date.today())) 
          historyPath = "%slogs_%s-%s/nose_library_build_history_%s" %(self.buildbotUrl, arch, system, compiler_version)
        elif str(test).find('test_example')!=-1 :
          (status,example) = list(test.test.arg)[:2]
          if (status == "build") : 
            path = example.path[example.path.rfind("/examples/")+1:]
          else :
            testPoint =  list(test.test.arg)[2]
            path = testPoint.path[testPoint.path.rfind("/examples/")+1:]
          arch = example.arch
          system = example.system
          compiler_version = example.compilerVersion
          logPath = "%slogs_%s-%s/%s/nose_%s_%s%s" %(self.buildbotUrl, arch, system, path, status, compiler_version, str(date.today())) 
          historyPath = "%slogs_%s-%s/%s/nose_%s_history_%s" %(self.buildbotUrl, arch, system, path, status, compiler_version)
        self.html.append('&nbsp;<a href="'+logPath+'">log</a>')
        self.html.append('&nbsp;<a href="'+historyPath+'">history</a>')
    
    def addSuccess(self, test):
        self.html.append('<a class="success">PASS</a>')
        self.insertLog(test)
        self.current=self.current.parent
        
    def addError(self, test, err):
        err = self.formatErr(err)
        self.html.append('<a class="fail">ERROR</a>')
        self.insertLog(test)
        self.current.setPass(False)
        self.current=self.current.parent
            
    def addFailure(self, test, err):
        err = self.formatErr(err)
        self.html.append('<a class="fail">FAIL</a>')
        self.insertLog(test)
        self.current.setPass(False)
        self.current=self.current.parent

    def finalize(self, result):
        for i in range(0,len(self.testLevelsInner)-1) :
          self.html.append('</ul>')
          if self.current.isPass():
            self.html.append('<a class="success">PASS</a>')                             
          else:
            self.html.append('<a class="fail">FAIL</a>')
          self.current=self.current.parent
          self.html.append('</li>')
        
        for i in range(0,len(self.testLevelsOuter)) :
          self.html.append('</ul>')
          if self.current.isPass():
            self.html.append('<a class="success">PASS</a>')                             
          else:
            self.html.append('<a class="fail">FAIL</a>')
          self.current=self.current.parent
          if i != len(self.testLevelsOuter)-1 :
            self.html.append('</li>')      
        if not result.wasSuccessful():
            self.html.append('<br /><a class="fail">'+str(len(result.failures))+' FAILs, '+str(len(result.errors))+ 'ERRORs</a>')                             
        else:
            self.html.append('<br /><a class="success">ALL PASS</a>')
        # print >> sys.stderr, self.html
        for l in self.html:
            self.stream.writeln(l)

    def formatErr(self, err):
        exctype, value, tb = err
        return ''.join(traceback.format_exception(exctype, value, tb))
    
    def setOutputStream(self, stream):
        # grab for own use
        self.stream = stream        
        return stream

    def startContext(self, ctx):
        try:
            n = ctx.__name__
            if n=='noseMain' and self.isStart :
              for l in self.header:
                self.stream.writeln(l)
              self.stream.writeln('<div style="display:none">')
              self.isStart=False
        except AttributeError:
            n = str(ctx).replace('<', '').replace('>', '')

    def stopContext(self, ctx):
        pass
    
    def startTest(self, test):
        description = ''
        ### OpenCMISS library build ###
        if str(test).find('test_build_library')!=-1 :
          description='Building the library'
          self.current = ResultTree(self.current)
        ### OpenCMISS example build, execute and check ###
        if str(test).find('test_example')!=-1 :
          (status, example) = list(test.test.arg)[:2]
          ## Build the example ##
          if status =='build' :
            for k in range(0,len(self.testLevelsInner)) :
              self.html.append('</ul>')
              if self.current.isPass():
                self.html.append('<a class="success">PASS</a>')                             
              else:
                self.html.append('<a class="fail">FAIL</a>')
              self.current=self.current.parent
              self.html.append('</li>')
            self.testLevelsInner = []
            path = example.path
            levels = path[path.rfind("/examples/")+1:].split('/')
            for i in range(0,len(levels)):
              if (len(self.testLevelsOuter)-1<i) :
                self.testLevelsOuter.append(levels[i])
                self.html.extend(['<li class="liClosed">&nbsp;',self.testLevelsOuter[i],'<ul>'])
                self.current = ResultTree(self.current)
              if (self.testLevelsOuter[i]!=levels[i]) :
                if(i<len(self.testLevelsOuter)-1) :
                  
                  for j in range(i,len(self.testLevelsOuter)-1) :
                    self.html.append('</ul>')
                    if self.current.isPass():
                      self.html.append('<a class="success">PASS</a>')                             
                    else:
                      self.html.append('<a class="fail">FAIL</a>')
                    self.current=self.current.parent
                    self.html.append('</li>')  
                self.testLevelsOuter[i]=levels[i]
                self.html.extend(['<li class="liClosed">&nbsp;',self.testLevelsOuter[i],'<ul>'])
                self.current = ResultTree(self.current)
                self.testLevelsOuter=self.testLevelsOuter[0:i+1]
            description='Building the test'
          ## Execute the example ##
          elif status=='run' :
            testPoint = list(test.test.arg)[2]
            path = testPoint.path
            if path!=example.path :
              levels = path[path.rfind(example.globalTestDir)+1:].split('/')
              for i in range(0,len(levels)):
                if (len(self.testLevelsInner)-1<i) :
                  self.testLevelsInner.append(levels[i])
                  self.html.extend(['<li class="liClosed">&nbsp;',self.testLevelsInner[i],'<ul>'])
                  self.current = ResultTree(self.current)
                if (self.testLevelsInner[i]!=levels[i]) :
                  if(i<len(self.testLevelsInner)-1) :
                    for j in range(i,len(self.testLevelsInner)-1) :
                      self.html.append('</ul>')
                      if self.current.isPass():
                        self.html.append('<a class="success">PASS</a>')                             
                      else:
                        self.html.append('<a class="fail">FAIL</a>')
                      self.current=self.current.parent
                      self.html.append('</li>')  
                  self.testLevelsInner[i]=levels[i]
                  self.html.extend(['<li class="liClosed">&nbsp;',self.testLevelsInner[i],'<ul>'])
                  self.current =  ResultTree(self.current)
                  self.testLevelsInner=self.testLevelsInner[0:i+1]
            description='Running the test'
          ## Check the output ##
          elif status=='check' :  
            description='Checking the output'
        self.current = ResultTree(self.current)
        self.html.extend([ '<li class="liBullet">&nbsp;',description,':&nbsp;'])
        
        
    def stopTest(self, test):
        self.html.append('</li>')
        if str(test).find('test_build_library')!=-1 :
          self.current=self.current.parent
        if str(test).find('test_example')!=-1 :
          status = list(test.test.arg)[0]
          if status=="run" or status=='check':
            testPoint = list(test.test.arg)[2]
            if (not hasattr(testPoint,'expectedPath')) or status=='check' :
              self.html.append('</ul>')
              if self.current.isPass():
                self.html.append('<a class="success">PASS</a>')                             
              else:
                self.html.append('<a class="fail">FAIL</a>')
              self.current=self.current.parent
              self.html.append('</li>')
         
import nose

if __name__ == '__main__':
    nose.main(addplugins=[HtmlOutput(this)])


