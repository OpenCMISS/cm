import urllib 
from BaseHTMLProcessor import BaseHTMLProcessor
                                      
sock = urllib.urlopen("../Programmer/html/examples.html") 
htmlSource = sock.read()                            
sock.close() 
parser = BaseHTMLProcessor()
parser.feed(htmlSource)    
f = open('../Programmer/html/examples1.html', 'w')
f.write(parser.output())                           