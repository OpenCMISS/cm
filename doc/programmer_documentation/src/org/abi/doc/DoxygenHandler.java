package org.abi.doc;

import java.io.IOException;
import java.io.StringWriter;
import java.util.Iterator;

import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import com.sun.org.apache.xml.internal.serialize.OutputFormat;
import com.sun.org.apache.xml.internal.serialize.XMLSerializer;

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
	private InterfaceElement interfaceElement;
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
			elementStarted = false;
			try {
				resultBuffer.append(getFunctionDocbookXML(functionElement));
			} catch (IOException e) {
				// TODO log
			}
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
				else if (name.equals("argsstring") )
				{
					functionElement.setArgs(textBuffer.toString().trim());
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
		if (name.equals("innerclass"))
		{
			interfaceElement.setName(textBuffer.toString().trim());
			elementStarted = false;
			try {
				resultBuffer.append(getInterfaceDocbookXML(interfaceElement));
			} catch (IOException e) {
				// TODO log
			}
			interfaceElement = null;
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
			if (name.equals("name") || name.equals("briefdescription") || name.equals("argsstring") 
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
		
		if (name.equals("innerclass"))
		{
			textBuffer = new StringBuffer();
			elementStarted = true;
			interfaceElement = new InterfaceElement();
		}
	}

	/**
	 * Write the xml results into docbook files. 
	 * @param element
	 * @return
	 */
	private String getFunctionDocbookXML(FunctionElement element) throws IOException, SAXException
	{
		OutputFormat of = new OutputFormat("XML","ISO-8859-1",true);
		of.setOmitXMLDeclaration(true);
		StringWriter writer = new StringWriter();
		XMLSerializer serializer = new XMLSerializer(writer,of);
		ContentHandler hd = serializer.asContentHandler();
		hd.startElement("","","command",null);
		addCharacters(hd, element.getName());
		addCharacters(hd, element.getArgs());
		hd.endElement("", "", "command");
		hd.startElement("", "", "itemizedlist", null);
		hd.startElement("", "", "listitem", null);
		addCharacters(hd, "Description: ");
		addCharacters(hd, element.getDescription());
		hd.endElement("", "", "listitem");
		hd.startElement("", "", "listitem", null);
		addCharacters(hd,"Parameters");
		hd.startElement("", "", "itemizedlist", null);
		Iterator<ParameterElement> iter = element.getParams().iterator();
		while(iter.hasNext())
		{
			ParameterElement pe = iter.next();
			hd.startElement("", "", "listitem", null);
			addCharacters(hd, pe.getName());
			addCharacters(hd, ": ");
			addCharacters(hd, pe.getDescription());
			hd.endElement("", "", "listitem");
		}
		hd.endElement("", "", "itemizedlist");
		hd.endElement("", "", "listitem");
		hd.endElement("", "", "itemizedlist");
		
		return writer.toString();
	}
	
	/**
	 * Write the xml results into docbook files. 
	 * @param element
	 * @return
	 */
	private String getInterfaceDocbookXML(InterfaceElement element) throws IOException, SAXException
	{
		OutputFormat of = new OutputFormat("XML","ISO-8859-1",true);
		of.setOmitXMLDeclaration(true);
		StringWriter writer = new StringWriter();
		XMLSerializer serializer = new XMLSerializer(writer,of);
		ContentHandler hd = serializer.asContentHandler();
		hd.startElement("","","command",null);
		hd.startElement("","","para",null);
		addCharacters(hd, element.getName());
		hd.endElement("", "", "para");
		hd.endElement("", "", "command");
		
		return writer.toString();
	}
	
	private void addCharacters(ContentHandler hd, String text) throws SAXException
	{
		if (text != null)
		{
			hd.characters(text.toCharArray(), 0, text.length());
		}
	}

}
