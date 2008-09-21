package org.abi.doc;

import java.util.ArrayList;
import java.util.List;

public class FunctionElement {
	private String name;
	private String description;
	private List<ParameterElement> params = new ArrayList<ParameterElement>();

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public boolean addParam(ParameterElement e) {
		return params.add(e);
	}

	public List<ParameterElement> getParams() {
		return params;
	}
}

class ParameterElement{
	private String type;
	private String description;
	private String name;
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
	}
	public String getDescription() {
		return description;
	}
	public void setDescription(String description) {
		this.description = description;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	
}
