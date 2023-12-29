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
 This class defines a contact interaction between a ligand atom and a residue atom.
 */

import java.io.*;
import java.util.*;

public class Contact {

	Residue res;
	/* Ligand atom*/
	Atom ref;
	/* Residue atom*/
	Atom trait;
	double distance;
	boolean selected;

	public Contact() {
	}

	public Contact(Residue res, Atom trait, Atom ref, double distance) {

		this.res = res;
		this.ref = ref;
		this.trait = trait;
		/* the distance is multiplied by 2 for better rendering in the diagram*/
		this.distance = 2 * distance;
		selected = false;

	}

	/* return the ligand atom*/
	public Atom getRef() {
		return ref;
	}

	/* return the residue atom*/
	public Atom getTrait() {
		return trait;
	}

	public Residue getResidue() {
		return res;
	}

	public double getDistance() {
		return distance;
	}

	public void setX(double x) {
		trait.setFx(x);
	}

	public void setY(double y) {
		trait.setFy(y);
	}

	public void setSelected(boolean b) {
		selected = b;
	}

	public boolean getSelected() {
		return selected;
	}

}
