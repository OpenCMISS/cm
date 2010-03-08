from buildbot.steps import shell
from buildbot.status.builder import SUCCESS, WARNINGS, FAILURE, SKIPPED, EXCEPTION

class ShellCommandWithHtmlLog(shell.ShellCommand):
     def evaluateCommand(self, cmd):
      if "failed" in self.getLog("results").getText():
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

class ShellCommandToCheckMissingRoutines(shell.ShellCommand):
     def evaluateCommand(self, cmd):
      if ("No functions missing in opencmiss_c.f90" in self.getLog("stdio").getText()) and ("No functions missing in opencmiss.h" in self.getLog("stdio").getText()):
        return SUCCESS
      else:
        return FAILURE

