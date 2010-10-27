"""This is a very basic example of a plugin that controls all test
output. In this case, it formats the output as ugly unstyled html.

Upgrading this plugin into one that uses a template and css to produce
nice-looking, easily-modifiable html output is left as an exercise for
the reader who would like to see his or her name in the nose AUTHORS file.
"""
import traceback
from nose.plugins import Plugin

import inspect

class HtmlOutput(Plugin):
    """Output test results as ugly, unstyled html.
    """
    
    name = 'html-output'
    score = 2 # run late
    runTestPass = True
    testPass = True
    all42TestPass = True
    testLevelsInner = []
    testLevelsOuter = []
    isFirst = True
    
    def __init__(self):
        super(HtmlOutput, self).__init__()
        self.html=['<html><head>',
                      '<title>Test output</title>',
                      '<SCRIPT LANGUAGE="JavaScript" SRC="tree.js"></SCRIPT>',
                      '<link href="tree.css" rel="stylesheet">',
                      '</head><body>',
                      '<ul id="tree1" class="mktree">']
    
    def addSuccess(self, test):
        self.html.append('<a class="success">PASS</a>')
        self.testPass = True
        if list(test.test.arg)[0]=='run'  :
          self.runTestPass = True   
        
    def addError(self, test, err):
        err = self.formatErr(err)
        self.html.append('<a class="fail">ERROR</a>')
        self.testPass = False
        self.all42TestPass = False
        if list(test.test.arg)[0]=='run'  :
          self.runTestPass = False 
            
    def addFailure(self, test, err):
        err = self.formatErr(err)
        self.html.append('<a class="fail">FAIL</a>')
        self.testPass = False
        self.all42TestPass = False
        if list(test.test.arg)[0]=='run'  :
          self.runTestPass = False

    def finalize(self, result):
        for i in range(0,len(self.testLevelsInner)-1) :
          self.html.append('</ul>')
          self.html.append('</li>')
        for i in range(0,len(self.testLevelsOuter)-1) :
          self.html.append('</ul>')
          if self.all42TestPass :
            self.html.append('<a class="success">PASS</a>')
          else :
            self.html.append('<a class="fail">FAIL</a>')
          self.html.append('</li>')
        self.html.append('</ul>')       
        if not result.wasSuccessful():
            self.html.append('<br><a class="fail">FAIL</a>')                             
        else:
            self.html.append('<br><a class="success">PASS</a>')
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
        except AttributeError:
            n = str(ctx).replace('<', '').replace('>', '')

    def stopContext(self, ctx):
        pass
    
    def startTest(self, test):
        description = ''
        if list(test.test.arg)[0]=='build' :
          path = list(test.test.arg)[1]
          path = path[path.find("/examples/")+10:]
          levels = path.split('/')
          for i in range(0,len(levels)):
            if (len(self.testLevelsOuter)-1<i) :
              self.testLevelsOuter.append(levels[i])
              self.html.extend(['<li class="liClosed">&nbsp;',self.testLevelsOuter[i],'<ul>'])
            if (self.testLevelsOuter[i]!=levels[i]) :
              if(i<len(self.testLevelsOuter)-1) :
                for k in range(0,len(self.testLevelsInner)) :
                  self.html.append('</ul>')
                  self.html.append('</li>')
                self.testLevelsInner = []
                for j in range(i,len(self.testLevelsOuter)-1) :
                  self.html.append('</ul>')
                  #TODO the checking of pass of fail is not correct.
                  if self.all42TestPass :
                    self.html.append('<a class="success">PASS</a>')
                  else :
                    self.html.append('<a class="fail">FAIL</a>')
                  self.html.append('</li>')  
                self.all42TestPass=True
              self.testLevelsOuter[i]=levels[i]
              self.html.extend(['<li class="liClosed">&nbsp;',self.testLevelsOuter[i],'<ul>'])
              self.testLevelsOuter=self.testLevelsOuter[0:i+1]
          description='Building the test'
        elif list(test.test.arg)[0]=='run' :
          path = list(test.test.arg)[6]
          levels = path.split('/')
          for i in range(0,len(levels)):
            if (len(self.testLevelsInner)-1<i) :
              self.testLevelsInner.append(levels[i])
              self.html.extend(['<li class="liClosed">&nbsp;',self.testLevelsInner[i],'<ul>'])
            if (self.testLevelsInner[i]!=levels[i]) :
              if(i<len(self.testLevelsInner)-1) :
                for j in range(i,len(self.testLevelsInner)-1) :
                  self.html.append('</ul>')
                  self.html.append('</li>')  
              self.testLevelsInner[i]=levels[i]
              self.html.extend(['<li class="liClosed">&nbsp;',self.testLevelsInner[i],'<ul>'])
              self.testLevelsInner=self.testLevelsInner[0:i+1]
          description='Running the test'
        elif list(test.test.arg)[0]=='check' :  
          description='Checking the output'
        self.html.extend([ '<li class="liBullet">&nbsp;',description,':&nbsp;'])
        
        
    def stopTest(self, test):
        self.html.append('</li>')
        if len(list(test.test.arg))==8 :
          self.html.append('</ul>')    
          if self.testPass:
            self.html.append('<a class="success">PASS</a>')                             
          else:
            self.html.append('<a class="fail">FAIL</a>')
          self.html.append('</li>')
        if list(test.test.arg)[0]=='check' :  
          self.html.append('</ul>')    
          if self.runTestPass and self.testPass:
            self.html.append('<a class="success">PASS</a>')                             
          else:
            self.html.append('<a class="fail">FAIL</a>')
          self.html.append('</li>')
         
import nose

if __name__ == '__main__':
    nose.main(addplugins=[HtmlOutput(this)])


