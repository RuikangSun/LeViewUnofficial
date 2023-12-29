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
 This class defines a water-mediated hydrogen bond.
 */

import java.util.*;

public class Path {

	Atom dep;
	Atom end;
	LinkedList waters;
	Atom coord;
	boolean selected;

	public Path() {
	}

	public Path(Atom a1, Atom a2, LinkedList l) {
		dep = a1;
		end = a2;
		waters = l;
		coord = a2.copy();
		selected = false;
	}

	public double getDist() {
		double r = 0;
		LinkedList tmp = new LinkedList();
		tmp.add(dep);
		for (int i = 0; i < waters.size(); i++) {
			Atom w = (Atom) waters.get(i);
			tmp.add(w);
		}
		tmp.add(end);

		for (int i = 0; i < tmp.size() - 1; i++) {
			Atom a1 = (Atom) tmp.get(i);
			Atom a2 = (Atom) tmp.get(i + 1);
			r = r + a1.distance(a2);
			i++;
		}

		return r;

	}

	public boolean getSelected() {
		return selected;
	}

	public void setSelected(boolean b) {
		selected = b;
	}

	public Atom getDep() {
		return dep;
	}

	public void setDep(Atom a) {
		dep = a;
	}

	public Atom getEnd() {
		return end;
	}

	public LinkedList getWaters() {
		return waters;
	}

	public Atom getCoord() {
		return coord;
	}

	public String getString() {

		String tmp = "";
		tmp = tmp + dep.getName() + " - ";

		for (int i = 0; i < waters.size(); i++) {
			Atom at = (Atom) waters.get(i);
			tmp = tmp + "H2O" + at.getId() + " ";
		}
		tmp = tmp + "- " + end.getName() + " (" + end.getParent().getName()
				+ end.getParent().getId() + ")";
		return tmp;

	}

	public String getString2() {

		String tmp = "";

		for (int i = 0; i < waters.size(); i++) {
			Atom at = (Atom) waters.get(i);
			tmp = tmp + "H2O" + at.getId() + ";";
		}
		tmp = tmp + end.getName() + " (" + end.getParent().getName()
				+ end.getParent().getId() + ")";
		return tmp;

	}

	public boolean isEqual(Path p) {
		if (p.getDep().getId() != dep.getId())
			return false;
		if (p.getEnd().getId() != end.getId())
			return false;

		LinkedList l = p.getWaters();

		if (l.size() != waters.size())
			return false;

		for (int i = 0; i < waters.size(); i++) {
			Atom a1 = (Atom) waters.get(i);
			Atom a2 = (Atom) l.get(i);
			if (a1.getId() != a2.getId())
				return false;
		}

		return true;

	}

	public void print() {
		System.out.println("le premier " + getDep() + "le second " + getEnd()
				+ "waters " + waters);
	}

}
