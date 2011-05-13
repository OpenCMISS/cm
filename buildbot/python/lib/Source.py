# Per-repository Scheduler

from buildbot.steps import source
from twisted.python import log


class SVNPropertyControlsMode(source.SVN):
    
    svn_revision = None    

    def __init__(self, **kwargs):
        if kwargs['workdir'].endswith('intel') :
            self.svn_revision = '36'
        source.SVN.__init__(self, **kwargs)

    def startVC(self, branch, revision, patch):
        try :
            self.args['mode'] = self.getProperty("SVNMode")
        except :
            print "SVNMode not defined, use the default mode"
        if self.args['mode'] == "clobber":
            self.description = ["checkout"]
            self.descriptionDone = ["checkout"]
        else:
            self.description = ["updating"]
            self.descriptionDone = ["update"]

        # Hack for opencmissextras to stick on older revision         
        if self.svn_revision is not None :
            revision = self.svn_revision 
            
        return source.SVN.startVC(self, branch, revision, patch)

#########################################
