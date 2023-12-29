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
 This class defines a ring.
 */

import java.io.*;
import java.util.*;

public class Ring {

	LinkedList atomId;
	boolean aromatic;
	LinkedList conjugate;
	boolean flag;
	double centerX;
	double centerY;
	double[] coord;
	
	public Ring() {
	}

	public Ring(LinkedList l) {
		atomId = l;
		aromatic = false;
		conjugate = new LinkedList();
		flag = false;
		centerX = 0;
		centerY = 0;
	}

	public void setAtoms(LinkedList ll) {
		atomId = ll;
	}

	public LinkedList getAtoms() {
		return atomId;
	}

	public double getCenterX() {
		return centerX;
	}

	public int getSize() {
		return atomId.size();
	}

	public double[] getCoord() {
		return coord;
	}

	public double getCenterY() {
		return centerY;
	}

	public void setCoord(double[] tab) {
		coord = tab;
	}

	public void setCenterX(double nx) {
		centerX = nx;
	}

	public void setCenterY(double ny) {
		centerY = ny;
	}

	public void setAromatic(boolean b) {
		aromatic = b;
	}

	public boolean getFlag() {
		return flag;
	}

	public void setFlag(boolean b) {
		flag = b;
	}

	public boolean isAromatic() {
		return aromatic;
	}

	public int getLenght() {
		return atomId.size();
	}

	public void addConjugate(Ring r) {
		conjugate.add(r);
	}

	public LinkedList getConjugate() {
		return conjugate;
	}

	public boolean isInside(LinkedList l) {
		for (int i = 0; i < atomId.size(); i++) {
			int id = (Integer) atomId.get(i);

			for (int j = 0; j < l.size(); j++) {
				int id2 = (Integer) l.get(j);
				if (id == id2)
					return true;
			}

		}
		return false;
	}

	public void print() {
		System.out.println("RING " + atomId + " " + aromatic + " " + conjugate);
	}

}
