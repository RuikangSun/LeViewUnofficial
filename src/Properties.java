/*
This file is part of LeView.

Copyright (C) 2013  Segolene Caboche

LeView is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeView is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeView.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 This class defines the properties used.
 */


import java.io.*;
import java.util.*;

public class Properties {

	LinkedList aa;
	LinkedList na;
	LinkedList ion;

	LinkedList elements;
	LinkedList radii;

	public Properties() {
	};

	public LinkedList getAA() {
		return aa;
	}

	public LinkedList getElements() {
		return elements;
	}

	public LinkedList getRadii() {
		return radii;
	}

	public LinkedList getNA() {
		return na;
	}

	public LinkedList getIon() {
		return ion;
	}

	public void setAA(LinkedList l) {
		this.aa = l;
	}

	public void setNA(LinkedList l) {
		this.na = l;
	}

	public void setElements(LinkedList l) {
		this.elements = l;
	}

	public void setRadii(LinkedList l) {
		this.radii = l;
	}

	public void setIon(LinkedList l) {
		this.ion = l;
	}

}
