package org.abi.doc;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;

import com.sun.org.apache.xml.internal.utils.DefaultErrorHandler;

/**
 * The main class to extract function information from doxygen xml files. parse it 
 * into relevant docbook xml format and insert it into docbook xml files.
 * @author tyu011
 */
public class DoxygenToDocbook {

	/**
	 * replace the relevant line in docbook xml file with info extracted from doxygen xml files.
	 * TODO Logging exceptions.
	 */
	public void replace(String doxygenPath, String doxygenFileName, StringBuffer result) {
		try {
			SAXParserFactory spf = SAXParserFactory.newInstance();
			spf.setNamespaceAware(true);
			SAXParser sp = spf.newSAXParser();
			InputSource input = new InputSource(new FileReader(new File(doxygenPath, doxygenFileName)));
			XMLReader reader = sp.getXMLReader(); 
			reader.setContentHandler(new DoxygenHandler(result));
			reader.setErrorHandler(new DefaultErrorHandler());
			reader.parse(input);
		} catch (ParserConfigurationException e) {
			// Log
		} catch (SAXException e) {
			// Log
		} catch (FileNotFoundException e) {
			// Log
		} catch (IOException e) {
			// Log
		}
	}

	/**
	 * TODO replace line identification be more flexible, e.g. distinguish with comments
	 * TODO use xml writer?
	 * @param args
	 */
	public static void main(String[] args) {
		String str;	
		try{
			FileInputStream	templateStream = new FileInputStream(args[0]);
			DataInputStream   input = new DataInputStream (templateStream);
			FileOutputStream outputStream = new FileOutputStream(args[0].substring(0, args[0].length()-9));
			DataOutputStream   output = new DataOutputStream(outputStream);

			while (null != ((str = input.readLine())))
			{
				int x=0;
				int y=0;
				StringBuffer result = new StringBuffer();
				while ((x=str.indexOf("<!--TO BE REPLACED ", y))>-1) {
					String moduleName = str.trim().substring(19, str.trim().length()-3);
					StringBuffer doxygenFileNameBuffer = new StringBuffer();
					doxygenFileNameBuffer.append("namespace");
					for(char a : moduleName.toCharArray())
					{
						doxygenFileNameBuffer.append('_');
						doxygenFileNameBuffer.append(a);
					}
					doxygenFileNameBuffer.append(".xml");
					result.append(str.substring(y,x));
					DoxygenToDocbook d2d = new DoxygenToDocbook();
					d2d.replace(args[1], doxygenFileNameBuffer.toString(), result);
					y = x + str.trim().length();
				}
				result.append(str.substring(y));
				str=result.toString();

				if(str.indexOf("'',") != -1){
					continue;
				}
				else{
					str=str+"\n";

					output.writeBytes(str);
				}
			}
		}
		catch (IOException ioe)
		{
			System.err.println ("I/O Error - " + ioe);
		}
	}
}

