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
import os

class ShellCommandWithHtmlTree(shell.ShellCommand):
     bcount = 0
     ecount = 0
     ccount = 0

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

     def operate(self,item,build,execute,check,parent=None):
       for i in list(item) : 
         if (i.tag=="li" and i.text.find("Building the test")!=-1) :
           if i.findtext("a")=="FAIL" :
             if (self.bcount%5==0 or self.btr==None) :
               self.btr=ET.Element("tr")
               build.append(self.btr)
             btd = ET.Element("td",width="20%")
             ba = list(i)[1]
             ba.text = parent.text
             btd.append(ba)
             self.btr.append(btd)
             self.bcount=self.bcount+1
         elif (i.tag=="li" and i.text.find("Running the test")!=-1) :
           if i.findtext("a")=="FAIL" :
             if (self.ecount%5==0 or self.etr==None) :
               self.etr=ET.Element("tr")
               execute.append(self.etr)
             etd = ET.Element("td",width="20%")
             ea = list(i)[1]
             ea.text = parent.text
             etd.append(ea)
             self.etr.append(etd)
             self.ecount=self.ecount+1
         elif (i.tag=="li" and i.text.find("Checking the output")!=-1) :
           if i.findtext("a")=="FAIL" :
             if (self.ccount%5==0 or self.ctr==None) :
               self.ctr=ET.Element("tr")
               check.append(self.ctr)
             ctd = ET.Element("td",width="20%")
             ca = list(i)[1]
             ca.text = parent.text
             ctd.append(ca)
             self.ctr.append(ctd)
             self.ccount=self.ccount+1
         else :
           self.operate(i,build,execute,check,item)


     def extractFails(self,content) :
       tree = ET.fromstring(content.replace("&nbsp;",""))
       ulelem = tree.find("body/ul")
       library_build = list(ulelem)[0]
       if (self.isFailed(library_build)) :
         outputContent = open("public_html/library_fail.html").read()
         outputtree = ET.fromstring(outputContent)
       else :
         outputContent = open("public_html/other_fails.html").read()
         outputtree = ET.fromstring(outputContent)
         for table in outputtree.findall("body/table") :
           if (table.get("id")=="execute") :
             outputexecute = table
           elif (table.get("id")=="build") :
             outputbuild = table
           elif (table.get("id")=="check") :
             outputcheck = table
         self.operate(ulelem, outputbuild,outputexecute,outputcheck)
         if self.bcount == 0 :
           outputtree.find("body").remove(outputbuild)
         if self.ecount == 0 :
           outputtree.find("body").remove(outputexecute)
         if self.ccount == 0 :
           outputtree.find("body").remove(outputcheck)
         self.addHTMLLog('fails', ET.tostring(outputtree))
   

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

