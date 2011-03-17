from buildbot.steps import shell
from buildbot.status.builder import SUCCESS, WARNINGS, FAILURE, SKIPPED, EXCEPTION

class ShellCommandWithHtmlLog(shell.ShellCommand):
     def evaluateCommand(self, cmd):
      if ("Fail" in self.getLog("results").getText()) :
        return FAILURE
      else:
        return SUCCESS

     def createSummary(self, log):
       stylesheet = """
       <style type="text/css">
       a.success {
       color: green;
       }
       a.fail {
       color: red;
       }
       </style>
       """
       myHeader = "<html><head>%s</head>" % (stylesheet,)
       myHeader += "<body>"
       content = myHeader
       content += log.getText()
       content += "</body></html>"

       self.addHTMLLog('results', content)

import xml.etree.ElementTree as ET
from xml.etree.ElementTree import XMLParser

class ShellCommandWithHtmlTree(shell.ShellCommand):
     def evaluateCommand(self, cmd):
      if ("FAIL" in self.getLog("results").getText()) :
        return FAILURE
      else:
        return SUCCESS

     def isFailed(self,item) :
       if item.text=="FAIL":
         return True
       else :
         for i in list(item) :
           if self.isFailed(i) :
             return True
         return False

     def operate(self,item):
       for i in list(item) : 
         if not self.isFailed(i) :
           if (not i.tag=="a") or i.text=="PASS" : 
             item.remove(i)
         else :
           self.operate(i)

     def extractFails(self,content) :
       tree = ET.fromstring(content.replace("&nbsp;",""))
       ulelem = tree.find("body/ul")
       self.operate(ulelem)
       self.addHTMLLog('fails', ET.tostring(tree))
   

     def createSummary(self, log):
       content =  log.getText()
       content += "</body></html>"
       self.addHTMLLog('results', content)
       self.extractFails(content)



class ShellCommandToCheckMissingRoutines(shell.ShellCommand):
     def evaluateCommand(self, cmd):
      if ("No functions missing in opencmiss_c.f90" in self.getLog("stdio").getText()) and ("No functions missing in opencmiss.h" in self.getLog("stdio").getText()):
        return SUCCESS
      else:
        return FAILURE

