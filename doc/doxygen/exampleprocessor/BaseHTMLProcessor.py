from sgmllib import SGMLParser
import htmlentitydefs

class BaseHTMLProcessor(SGMLParser):

	def reset(self):
		# extend (called by SGMLParser.__init__)
		self.pieces = []
		self.buffer = ""
		self.parents = []
		SGMLParser.reset(self)
		
	def unknown_starttag(self, tag, attrs):
		strattrs = "".join([' %s="%s"' % (key, value) for key, value in attrs])
		if tag=='a' :
			self.buffer = "<%(tag)s%(strattrs)s>" % locals()
		elif tag != 'li':
			self.pieces.append("<%(tag)s%(strattrs)s>" % locals())
		
	def unknown_endtag(self, tag):
		if tag != 'li' :
			self.pieces.append("</%(tag)s>" % locals())

	def handle_charref(self, ref):
		self.pieces.append("&#%(ref)s;" % locals())
		
	def handle_entityref(self, ref):
		self.pieces.append("&%(ref)s" % locals())
		if htmlentitydefs.entitydefs.has_key(ref):
			self.pieces.append(";")

	def handle_data(self, text):
		rmindex=text.find('/src/')
		if rmindex>0:
			text=text[:rmindex]
			texts = text.split('/')
			for i in range (0,len(self.parents)) :
			  if(texts[i]!=self.parents[i]) :
			    for j in range (i, len(self.parents)) :
			  		self.pieces.append("</li></ul>")
			    self.parents = self.parents[:i]
			    break
			self.pieces.append("<li>")
			for i in range (len(self.parents),len(texts)-1) :
				self.parents.append(texts[i])
				self.pieces.append(texts[i]+"<ul><li>")
			self.pieces.append(self.buffer)
			self.buffer=""
			self.pieces.append(texts[len(texts)-1])
		else :
			self.pieces.append(text)
		
	def handle_comment(self, text):
		self.pieces.append("<!--%(text)s-->" % locals())
		
	def handle_pi(self, text):
		self.pieces.append("<?%(text)s>" % locals())

	def handle_decl(self, text):
		self.pieces.append("<!%(text)s>" % locals())
		
	def output(self):
		"""Return processed HTML as a single string"""
		return "".join(self.pieces)

if __name__ == "__main__":
	for k, v in globals().items():
		print k, "=", v
