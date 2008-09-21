package org.abi.doc;

import java.util.Iterator;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Read doxygen xml files using SAX. Extract functions out.
 * TODO allow to extract interface as well.
 * 
 * @author tyu011
 */
public class DoxygenHandler extends DefaultHandler {

	private boolean functionStarted = false;
	private boolean elementStarted = false;
	private boolean paramStarted = false;

	private FunctionElement functionElement;
	private ParameterElement parameterElement;

	private StringBuffer textBuffer;
	private StringBuffer resultBuffer;


	public DoxygenHandler(StringBuffer result) {
		resultBuffer = result;
	}

	@Override
	public void characters(char[] ch, int start, int length)
	throws SAXException {
		if (elementStarted)
		{
			textBuffer.append(ch, start, length);
		}
	}

	@Override
	public void endElement(String uri, String localName, String name)
	throws SAXException {
		if (name.equals("memberdef") && functionStarted)
		{
			functionStarted = false;
			resultBuffer.append(getFunctionDocbookXML(functionElement));
			functionElement = null;
		}
		if (functionStarted)
		{
			if (!paramStarted)
			{
				if (name.equals("name"))
				{
					functionElement.setName(textBuffer.toString().trim());
					elementStarted = false;
				}
				else if (name.equals("briefdescription") )
				{
					functionElement.setDescription(textBuffer.toString().trim());
					elementStarted = false;
				}
			}
			else
			{
				if (name.equals("param"))
				{
					paramStarted = false;
					if (parameterElement.getName() == null)
					{
						parameterElement.setName("*");
					}
					functionElement.addParam(parameterElement);
				}
				else if (name.equals("defname"))
				{
					parameterElement.setName(textBuffer.toString().trim());
					elementStarted = false;
				}
				else if (name.equals("briefdescription") )
				{
					parameterElement.setDescription(textBuffer.toString().trim());
					elementStarted = false;
				}
				else if (name.equals("type") )
				{
					parameterElement.setType(textBuffer.toString().trim());
					elementStarted = false;
				}	
			}

		}

	}

	@Override
	public void startElement(String uri, String localName, String name,
			Attributes atts) throws SAXException {
		if (name.equals("memberdef") && "function".equals(atts.getValue("kind")) 
				&& "public".equals(atts.getValue("prot")))
		{
			functionStarted = true;
			functionElement = new FunctionElement();
		}
		if (functionStarted)
		{
			// Inside a function
			if (name.equals("name") || name.equals("briefdescription") 
					|| name.equals("type") || name.equals("defname"))
			{
				textBuffer = new StringBuffer();
				elementStarted = true;
			}
			else if (name.equals("param"))
			{
				paramStarted = true;
				parameterElement = new ParameterElement();
			}
		}
	}

	/**
	 * Write the xml results into docbook files. 
	 * TODO use xml writer?
	 * @param element
	 * @return
	 */
	private String getFunctionDocbookXML(FunctionElement element)
	{
		StringBuffer buf = new StringBuffer();
		buf.append("      <command>");
		buf.append(element.getName());
		buf.append("</command>\n");
		buf.append("      <itemizedlist>\n");
		buf.append("        <listitem>Description:");
		buf.append(element.getDescription() != null ? element.getDescription() : "");
		buf.append("</listitem>\n");
		buf.append("        <listitem>Parameters:\n");
		buf.append("          <itemizedlist>\n");
		Iterator<ParameterElement> iter = element.getParams().iterator();
		while(iter.hasNext())
		{
			ParameterElement pe = iter.next();
			buf.append("            <listitem>");

			buf.append(pe.getName());
			buf.append(": ");
			buf.append(pe.getDescription()!=null ? escapeXML(pe.getDescription()) : "");
			buf.append("</listitem>\n");
		}
		buf.append("          </itemizedlist>\n");
		buf.append("        </listitem>\n");
		buf.append("      </itemizedlist>\n");
		return buf.toString();
	}

	private String escapeXML(String text)
	{
		String newText = text.replace("<", "&lt;");
		newText = newText.replace(">", "&gt;");
		return newText;
	}

}
