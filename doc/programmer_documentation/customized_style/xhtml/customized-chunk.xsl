<?xml version="1.0" encoding="US-ASCII"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" xmlns="http://www.w3.org/1999/xhtml" version="1.0" exclude-result-prefixes="exsl">

<xsl:import href="../../nwalsh/xhtml/chunk.xsl"/>

<xsl:param name="html.ext" select="'.xhtml'"/>
  
<xsl:template name="section.titlepage.before.recto">
  <!-- TODO This part is the customised part for the return to top link -->
  <xsl:variable name="level">
    <xsl:call-template name="section.level"/>
  </xsl:variable>
  <xsl:variable name="chunkfn">
    <xsl:apply-templates mode="chunk-filename" select="."/>
  </xsl:variable>

  <xsl:if test="$level &gt; $chunk.section.depth">
    <p class="returntotop">
      <a href="{$chunkfn}">
        <xsl:text>Return to top</xsl:text>
      </a>
    </p>  
  </xsl:if>
  <!-- END of customised part-->
</xsl:template>  


</xsl:stylesheet>
