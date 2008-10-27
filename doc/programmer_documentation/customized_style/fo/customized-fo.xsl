<?xml version="1.0" encoding="US-ASCII"?>
<xsl:stylesheet 
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
  xmlns:exsl="http://exslt.org/common" 
  xmlns="http://www.w3.org/1999/xhtml" 
  xmlns:fo="http://www.w3.org/1999/XSL/Format"
  xmlns:mml="http://www.w3.org/1998/Math/MathML" version="1.0" 
  exclude-result-prefixes="exsl">

<xsl:import href="../../nwalsh/fo/docbook.xsl"/>
  
<xsl:param name="body.start.indent" select="'0pt'"/>
<xsl:attribute-set name="section.title.level1.properties">
  <xsl:attribute name="font-size">
    <xsl:value-of select="$body.font.master * 1.4"/>
    <xsl:text>pt</xsl:text>
  </xsl:attribute>
</xsl:attribute-set>
<xsl:attribute-set name="section.title.level2.properties">
  <xsl:attribute name="font-size">
    <xsl:value-of select="$body.font.master * 1.3"/>
    <xsl:text>pt</xsl:text>
  </xsl:attribute>
</xsl:attribute-set>
<xsl:attribute-set name="section.title.level3.properties">
  <xsl:attribute name="font-size">
    <xsl:value-of select="$body.font.master * 1.2"/>
    <xsl:text>pt</xsl:text>
  </xsl:attribute>
</xsl:attribute-set>
<xsl:attribute-set name="section.title.level4.properties">
  <xsl:attribute name="font-size">
    <xsl:value-of select="$body.font.master * 1.1"/>
    <xsl:text>pt</xsl:text>
  </xsl:attribute>
</xsl:attribute-set>
  
  
<xsl:template match="inlineequation/mml:math">
  <fo:instream-foreign-object>
  <xsl:copy>
    <xsl:apply-templates mode="copy-all"/>
  </xsl:copy>
  </fo:instream-foreign-object>
</xsl:template>  


</xsl:stylesheet>
