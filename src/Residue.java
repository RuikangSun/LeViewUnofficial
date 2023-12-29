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
 This class defines a residue.
 */

import java.io.*;
import java.util.*;

public class Residue {

	String name;
	LinkedList atoms;
	String id;
	Chain parent;
	String type;
	LinkedList neig;
	LinkedList bonds;
	LinkedList Hneig;
	boolean mapping;

	double[] Mapping = new double[4];

	public Residue() {
		atoms = new LinkedList();
		neig = new LinkedList();
		bonds = new LinkedList();
		Hneig = new LinkedList();
		mapping = false;
	}

	public String getName() {
		return name;
	}

	public Residue copy() {
		Residue r = new Residue();
		r.setName(name);
		r.setId(id);
		r.setType(type);
		r.setParent(parent);
		r.setHneig(Hneig);
		r.setNeig(neig);

		LinkedList l = new LinkedList();
		for (int i = 0; i < atoms.size(); i++) {
			Atom a = (Atom) atoms.get(i);
			l.add(a.copy());
		}
		r.setAtoms(l);

		return r;

	}

	public LinkedList getAtoms() {
		return atoms;
	}

	public String getId() {
		return id;
	}

	public void setMapping(boolean b) {
		mapping = b;
	}

	public boolean getMapping() {
		return mapping;
	}

	public void addResNeig(String ID) {
		neig.add(ID);
	}

	public void addResHneig(String ID) {
		Hneig.add(ID);
	}

	public LinkedList getNeig() {
		return neig;
	}

	public LinkedList getHneig() {
		return Hneig;
	}

	public void setHneig(LinkedList l) {
		Hneig = l;
	}

	public void setNeig(LinkedList l) {
		neig = neig;
	}

	public int getIdInt() {
		String tp1 = getId();
		String ttp1 = "";
		for (int i = 0; i < tp1.length(); i++) {
			char cc = tp1.charAt(i);
			if (Character.isDigit(cc)) {
				ttp1 = ttp1 + cc;
			}
		}
		return Integer.parseInt(ttp1);
	}

	public void setName(String c) {
		name = c;
	}

	public Chain getParent() {
		return parent;
	}

	public void setId(String i) {
		id = i;
	}

	public void setParent(Chain c) {
		parent = c;
	}

	public void addAtom(Atom a) {
		atoms.add(a);
	}

	public void setAtoms(LinkedList l) {
		atoms = l;
	}

	public void setType() {
		type = type();
	}

	public void setType(String tt) {
		type = tt;
	}

	public String getType() {
		return type;
	}

	public boolean isIon2() {
		if ((name.equals("MG")) || (name.equals("FE")) || (name.equals("ZN"))
				|| (name.equals("CA")) || (name.equals("HG"))
				|| (name.equals("MN")) || (name.equals("CL"))
				|| (name.equals("BR")) || (name.equals("I"))
				|| (name.equals("IOD")) || (name.equals("CD"))
				|| (name.equals("K")) || (name.equals("NA"))
				|| (name.equals("NI")) || (name.equals("CU")))
			return true;
		else
			return false;
	}

	public boolean isIon(Chain chain) {
		Properties prop = chain.getParent().getProp();
		LinkedList ion = prop.getIon();
		if (ion.indexOf(name) != -1)
			return true;
		return false;

	}

	public double getDistMax(double x, double y) {

		double dist = 0;

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);

			double d = Math.sqrt((x - at.getFx()) * (x - at.getFx())
					+ (y - at.getFy()) * (y - at.getFy()));
			if (d > dist) {
				dist = d;
			}
		}

		return dist;

	}

	public double[] getFront() {

		double[] tab = new double[4];

		double xmin = 10000;
		double xmax = -10000;
		double ymin = 10000;
		double ymax = -10000;

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);

			if (at.getFx() < xmin)
				xmin = at.getFx();
			if (at.getFx() > xmax)
				xmax = at.getFx();
			if (at.getFy() < ymin)
				ymin = at.getFx();
			if (at.getFy() < ymax)
				ymax = at.getFx();

		}

		tab[0] = xmin;
		tab[1] = xmax;
		tab[2] = ymin;
		tab[3] = ymax;

		return tab;

	}

	public int getNatom() {
		return atoms.size();
	}

	public int getHetero() {
		int c = 0;
		for (int i = 0; i < atoms.size(); i++) {
			Atom a = (Atom) atoms.get(i);
			if (a.getHetero())
				c++;
		}
		return c;
	}

	public boolean isLinked(Residue r) {
		LinkedList atoms2 = r.getAtoms();
		for (int i = 0; i < atoms.size(); i++) {
			Atom at1 = (Atom) atoms.get(i);
			for (int j = 0; j < atoms2.size(); j++) {
				Atom at2 = (Atom) atoms2.get(j);
				if (at1.isLinked(at2))
					return true;
			}
		}
		return false;
	}

	public String type2() {
		if ((name.equals("ALA")) || (name.equals("ARG"))
				|| (name.equals("ASN")) || (name.equals("ASP"))
				|| (name.equals("CYS")) || (name.equals("GLU"))
				|| (name.equals("GLN")) || (name.equals("GLY"))
				|| (name.equals("HIS")) || (name.equals("ILE"))
				|| (name.equals("LEU")) || (name.equals("LYS"))
				|| (name.equals("MET")) || (name.equals("PHE"))
				|| (name.equals("PRO")) || (name.equals("SER"))
				|| (name.equals("TRP")) || (name.equals("THR"))
				|| (name.equals("TYR")) || (name.equals("VAL")))
			return "aa";
		if ((name.equals("GUA")) || (name.equals("CYT"))
				|| (name.equals("ADE")) || (name.equals("THY"))
				|| (name.equals("G")) || (name.equals("U"))
				|| (name.equals("A")) || (name.equals("C"))
				|| (name.equals("T")) || (name.equals("DA"))
				|| (name.equals("DC")) || (name.equals("DG"))
				|| (name.equals("DT")) || (name.equals("5CM"))
				|| (name.equals("NRI")) || (name.equals("8OG")))
			return "dna";
		if ((name.equals("MG")) || (name.equals("FE")) || (name.equals("ZN"))
				|| (name.equals("CA")) || (name.equals("HG"))
				|| (name.equals("MN")) || (name.equals("CL"))
				|| (name.equals("BR")) || (name.equals("I"))
				|| (name.equals("IOD")) || (name.equals("CD"))
				|| (name.equals("K")) || (name.equals("NA"))
				|| (name.equals("NI")) || (name.equals("CU")))
			return "ion";
		for (int i = 0; i < atoms.size(); i++) {
			Atom a = (Atom) atoms.get(i);
			if (a.getName().equals("CA"))
				return "ns";
		}
		return "unknown";
	}

	public String type() {
		Properties prop = parent.getParent().getProp();
		LinkedList aa = prop.getAA();
		LinkedList na = prop.getNA();
		LinkedList ion = prop.getIon();
		if (na.indexOf(name) != -1)
			return "dna";
		if (aa.indexOf(name) != -1)
			return "aa";
		if (ion.indexOf(name) != -1)
			return "ion";
		for (int i = 0; i < atoms.size(); i++) {
			Atom a = (Atom) atoms.get(i);
			if (a.getName().equals("CA"))
				return "ns";
		}
		return "unknown";
	}

	public Atom getAtom(int ID) {
		for (int i = 0; i < atoms.size(); i++) {
			Atom a = (Atom) atoms.get(i);
			if (a.getId() == ID)
				return a;
		}
		return null;
	}

	public int getAtomPos(int ID) {
		for (int i = 0; i < atoms.size(); i++) {
			Atom a = (Atom) atoms.get(i);
			if (a.getId() == ID)
				return i;
		}
		return -1;
	}

	public void print() {
		System.out.println("########################");
		System.out.println("RESIDUE " + id + " " + name + " " + type + " "
				+ neig + " " + mapping);
		System.out.println("number of atoms " + atoms.size());
		for (int i = 0; i < atoms.size(); i++) {
			Atom c = (Atom) atoms.get(i);
			c.print();
		}
		System.out.println("########################");
	}

	public final void setMapping(double[] map) {
		Mapping = map;
	}

}
