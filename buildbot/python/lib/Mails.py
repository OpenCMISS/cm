from buildbot.status.mail import MailNotifier
from buildbot.status.builder import Results
from email.Utils import formatdate
from email.MIMEText import MIMEText
from email.MIMEMultipart import MIMEMultipart

ENCODING = 'utf8'

class MailNotifierWithHtmlAttachment(MailNotifier):
     def createEmail(self, msgdict, builderName, projectName, results, 
                    patch=None, logs=None):
        text = msgdict['body'].encode(ENCODING)
        type = msgdict['type']
        if 'subject' in msgdict:
            subject = msgdict['subject'].encode(ENCODING)
        else:
            subject = self.subject % { 'result': Results[results],
                                       'projectName': projectName,
                                       'builder': builderName,
                                       }


        assert type in ('plain', 'html'), "'%s' message type must be 'plain' or 'html'." % type

        if patch or logs:
            m = MIMEMultipart()
            m.attach(MIMEText(text, type, ENCODING))
        else:
            m = Message()
            m.set_payload(text, ENCODING)
            m.set_type("text/%s" % type)

        m['Date'] = formatdate(localtime=True)
        m['Subject'] = subject
        m['From'] = self.fromaddr
        # m['To'] is added later

        if patch:
            a = MIMEText(patch[1].encode(ENCODING), _charset=ENCODING)
            a.add_header('Content-Disposition', "attachment",
                         filename="source patch")
            m.attach(a)
        if logs:
            for log in logs:
                name = "%s.%s" % (log.getStep().getName(),
                                  log.getName())
                if self._shouldAttachLog(log.getName()) or self._shouldAttachLog(name):
                    a = MIMEText(log.getText().encode(ENCODING),_subtype="html",
                                 _charset=ENCODING)
                    a.add_header('Content-Disposition', "attachment",
                                 filename=name)
                    m.attach(a)

        # Add any extra headers that were requested, doing WithProperties
        # interpolation if necessary
        if self.extraHeaders:
            for k,v in self.extraHeaders.items():
                k = properties.render(k)
                if k in m:
                    twlog("Warning: Got header " + k + " in self.extraHeaders "
                          "but it already exists in the Message - "
                          "not adding it.")
                    continue
                m[k] = properties.render(v)

        return m
