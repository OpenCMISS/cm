<?xml version="1.0" encoding="US-ASCII"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" xmlns="http://www.w3.org/1999/xhtml" version="1.0" exclude-result-prefixes="exsl">

<xsl:import href="../../nwalsh/xhtml/docbook.xsl"/>

<xsl:param name="html.ext" select="'.xhtml'"/>
  
<xsl:template name="section.titlepage.before.recto">
  <xsl:variable name="top-anchor">
    <xsl:call-template name="object.id">
      <xsl:with-param name="object" select="/*[1]"/>
    </xsl:call-template>
  </xsl:variable>

    <p class="returntotop">
      <a href="#{$top-anchor}">
        <xsl:text>Return to top</xsl:text>
      </a>
    </p>
</xsl:template>
  
</xsl:stylesheet>