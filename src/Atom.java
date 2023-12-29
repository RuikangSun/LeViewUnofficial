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
This class defines an atom.
 */

import java.io.*;
import java.util.*;

public class Atom {

	int waterId;
	int type;
	/*
	 * 1 chaine 2 pair 3 ring 4 aromatic ring 5 isolated 6 terminal
	 */
	String name;
	int id;
	String element;
	double x;
	double y;
	double z;
	Residue parent;
	boolean hetero;
	/* List of connected atoms */
	LinkedList neig;
	/* List of atom having an explicit bond with this one */
	LinkedList Hneig;
	int wType;
	LinkedList wMediated;
	/* Number of hydrogens */
	int Hnumber;
	/* true if this atom is a terminal one */
	boolean terminal;
	/* 2D coordinates */
	double fx;
	double fy;
	/* true if this has 2D coordinates */
	boolean flag;
	/*
	 * true if this atom can be use as a pivot atom for resolving conflicts in
	 * ligand layout
	 */
	boolean pivot;
	boolean pivot2;
	int xg;
	int yg;
	int xoff = 0;
	int yoff = 0;
	/* true if this atom is selected in the graphical interface */
	boolean selected;

	public Atom() {
		/* initialisation */
		neig = new LinkedList();
		Hneig = new LinkedList();
		wMediated = new LinkedList();
		LinkedList l1 = new LinkedList();
		wMediated.add(l1);
		LinkedList l2 = new LinkedList();
		wMediated.add(l2);
		LinkedList l3 = new LinkedList();
		wMediated.add(l3);
		Hnumber = 0;
		terminal = false;
		type = -1;
		wType = -1;
		flag = false;
		pivot = false;
		pivot2 = false;
		int xg = 0;
		int yg = 0;
		int waterId = 0;
		selected = false;
	}

	public void setTerm(boolean b) {
		terminal = b;
	}

	public void setSelected(boolean b) {
		selected = b;
	}

	public boolean getSelected() {
		return selected;
	}

	/* initailisation of some parameters */
	public void init() {
		setFx(0);
		setFy(0);
		Hnumber = 0;
		terminal = false;
		type = -1;
		flag = false;
		wMediated = new LinkedList();
		LinkedList l1 = new LinkedList();
		wMediated.add(l1);
		LinkedList l2 = new LinkedList();
		wMediated.add(l2);
		LinkedList l3 = new LinkedList();
		wMediated.add(l3);
	}

	public void setWaterId(int i) {
		waterId = i;
	}

	public int getWaterId() {
		return waterId;
	}

	public void setXg(int i) {
		xg = i;
	}

	public void setYg(int i) {
		yg = i;
	}

	public int getXg() {
		return xg;
	}

	public int getYg() {
		return yg;
	}

	public int getwType() {
		return wType;
	}

	public void setwType(int i) {
		wType = i;
	}

	public boolean getPivot() {
		return pivot;
	}

	public boolean getPivot2() {
		return pivot2;
	}

	public void setPivot(boolean b) {
		pivot = b;
	}

	public void setPivot2(boolean b) {
		pivot2 = b;
	}

	public LinkedList getwMediated() {
		return wMediated;
	}

	public LinkedList getwMediated1() {
		LinkedList l = (LinkedList) wMediated.get(0);
		return l;
	}

	public LinkedList getwMediated2() {
		LinkedList l = (LinkedList) wMediated.get(1);
		return l;
	}

	public LinkedList getwMediated3() {
		LinkedList l = (LinkedList) wMediated.get(2);
		return l;
	}

	public void addwMediated1(Atom a) {
		LinkedList l = (LinkedList) wMediated.get(0);
		l.add(a);
	}

	public void addwMediated2(Atom a) {
		LinkedList l = (LinkedList) wMediated.get(1);
		l.add(a);
	}

	public void addwMediated3(Atom a) {
		LinkedList l = (LinkedList) wMediated.get(2);
		l.add(a);
	}

	public Atom clone() {
		Atom a = new Atom();
		a.setNeig(neig);
		a.setHneig(Hneig);
		a.setFlag(flag);
		a.setName(name);
		a.setFx(fx);
		a.setFy(fy);
		a.setTerm(terminal);
		return a;
	}

	public Atom copy() {
		Atom a = new Atom();
		a.setNeig(neig);
		a.setHneig(Hneig);
		a.setName(name);
		a.setElement(element);
		a.setTerm(terminal);
		a.setParent(parent);
		a.setCoord(x, y, z);
		a.setId(id);
		return a;
	}

	public Atom getCopy() {

		Atom a = new Atom();
		a.setNeig(neig);
		a.setHneig(Hneig);
		a.setFlag(flag);
		a.setName(name);
		a.setFx(fx);
		a.setFy(fy);
		a.setElement(element);
		a.setTerm(terminal);
		a.setParent(parent);
		a.setCoord(x, y, z);
		a.setId(id);
		return a;

	}

	public boolean getFlag() {
		return flag;
	}

	public void setFlag(boolean b) {
		flag = b;
	}

	public boolean getTerm() {
		return terminal;
	}

	public void setType(int type) {
		this.type = type;
	}

	public int getType() {
		return type;
	}

	public String getName() {
		return name;
	}

	public void setFx(double d) {
		fx = d;
	}

	public void setFy(double d) {
		fy = d;
	}

	public double getFx() {
		return fx;
	}

	public double getFy() {
		return fy;
	}

	public String getElement() {
		return element;
	}

	public int getId() {
		return id;
	}

	public int getHnumber() {
		return Hnumber;
	}

	public void setHnumber(int h) {
		Hnumber = h;
	}

	public Residue getParent() {
		return parent;
	}

	public double getX() {
		return x;
	}

	public double getY() {
		return y;
	}

	public double getZ() {
		return z;
	}

	/* Used for scale function in the graphical interface */
	public final void getSize(double[] mapping, Atom[] view) {
		if (fx < mapping[0]) {
			mapping[0] = fx;
			view[0] = this;
		}
		if (fx > mapping[1]) {
			mapping[1] = fx;
			view[1] = this;
		}
		if (fy < mapping[2]) {
			mapping[2] = fy;
			view[2] = this;
		}
		if (fy > mapping[3]) {
			mapping[3] = fy;
			view[3] = this;
		}
	}

	public final int getFX(double[] mapping, double bit) {
		return (int) (mapping[0] * (fx - mapping[1]) + bit) + xoff;
	}

	public final int getFY(double[] mapping, double bit) {
		return (int) (mapping[0] * (fy - mapping[2]) + bit) + yoff;
	}

	public void setName(String s) {
		name = s;
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

	public int Nneig() {
		return neig.size();
	}

	public void setElement(String s) {
		element = s;
	}

	public void setId(int i) {
		id = i;
	}

	public void setParent(Residue r) {
		parent = r;
	}

	public void setCoord(double a, double b, double c) {
		x = a;
		y = b;
		z = c;
	}

	public void setHetero(boolean f) {
		hetero = f;
	}

	public boolean getHetero() {
		return hetero;
	}

	public void setNeig(LinkedList l) {
		neig = l;
	}

	public boolean addNeig(Atom ID) {
		for (int i = 0; i < neig.size(); i++) {
			Atom at = (Atom) neig.get(i);
			if (at.getId() == ID.getId())
				return false;
		}
		neig.add(ID);
		return true;
	}

	public boolean addHneig(Atom ID) {
		for (int i = 0; i < Hneig.size(); i++) {
			Atom at = (Atom) Hneig.get(i);
			if (at.getId() == ID.getId())
				return false;
		}

		for (int i = 0; i < neig.size(); i++) {
			Atom at = (Atom) neig.get(i);
			if (at.getId() == ID.getId())
				return false;
		}

		Hneig.add(ID);
		return true;
	}

	/* return true if the parameter atom has a covalent bond with this one */
	public boolean isLinked(Atom at) {
		if (id == at.getId())
			return false;
		for (int i = 0; i < neig.size(); i++) {
			Atom n = (Atom) neig.get(i);
			if (n.getId() == at.getId())
				return true;
		}
		return false;
	}

	/* return the distance from 3D coordinates */
	public double distance(Atom at) {
		double dx = x - at.getX();
		double dy = y - at.getY();
		double dz = z - at.getZ();

		double dd = dx * dx + dy * dy + dz * dz;
		if (dd > 0)
			dd = Math.sqrt(dd);
		else
			dd = 0.0F;
		return dd;
	}

	/* return the distance from the 2D coordinates */
	public double distanceF(Atom at) {
		double dx = fx - at.getFx();
		double dy = fy - at.getFy();

		double dd = dx * dx + dy * dy;
		if (dd > 0)
			dd = Math.sqrt(dd);
		else
			dd = 0.0F;
		return dd;
	}

	public void removeNeig(int e) {
		int c = 0;
		for (int i = 0; i < neig.size(); i++) {
			Atom att = (Atom) neig.get(i);
			if (att.getId() == e) {
				break;
			} else {
				c++;
			}
		}
		neig.remove(c);
	}

	/**
	 * return true if this atom is an H-Bond donor
	 */
	public boolean isHBondDonor() {
		if (!parent.getType().equals("aa")) {
			return false;
		} else {
			String resname = parent.getName();

			if (name.equals("N")) {
				if ("PRO".equals(resname)) {
					return false;
				} else {
					return true;
				}
			} else if (name.equals("NE1") && "TRP".equals(resname)) {
				return true;
			} else if (name.equals("OG") && "SER".equals(resname)) {
				return true;
			} else if (name.equals("OG1") && "THR".equals(resname)) {
				return true;
			} else if (name.equals("ND2") && "ASN".equals(resname)) {
				return true;
			} else if (name.equals("NE2") && "GLN".equals(resname)) {
				return true;
			} else if (name.equals("OH") && "TYR".equals(resname)) {
				return true;
			} else if (name.equals("NE1") && "HIS".equals(resname)) {
				return true;
			} else if ((name.equals("NE2") || name.equals("ND1"))
					&& "HIS".equals(resname)) {
				return true;
			} else if (name.equals("NZ") && "LYS".equals(resname)) {
				return true;
			} else if ((name.equals("NE") || name.equals("NH1") || name
					.equals("NH2"))
					&& "ARG".equals(resname)) {
				return true;
			} else {
				return false;
			}
		}
	}

	public boolean isHBondDonor2() {
		if ((element.equals("O")) && (Hnumber != 0))
			return true;
		if ((element.equals("N")) && (Hnumber != 0))
			return true;
		if ((element.equals("S")) && (Hnumber != 0))
			return true;
		return false;
	}

	public boolean isHBondAcceptor2() {
		if (element.equals("O"))
			return true;
		if (element.equals("S"))
			return true;
		if (element.equals("N"))
			return true;
		return false;
	}

	/**
	 * return true if this atom is an H-Bond acceptor
	 */
	public boolean isHBondAcceptor() {

		if (!parent.getType().equals("aa")) {
			return false;
		} else if (name.equals("O")) {
			return true;
		} else {
			String resname = parent.getName();

			if (name.equals("NE1") && "TRP".equals(resname)) {
				return true;
			} else if (name.equals("OG") && "SER".equals(resname)) {
				return true;
			} else if (name.equals("OG1") && "THR".equals(resname)) {
				return true;
			} else if (name.equals("OD1") && "ASN".equals(resname)) {
				return true;
			} else if (name.equals("OE1") && "GLN".equals(resname)) {
				return true;
			} else if (name.equals("OH") && "TYR".equals(resname)) {
				return true;
			} else if (name.equals("NE1") && "HIS".equals(resname)) {
				return true;
			} else if ((name.equals("NE2") || name.equals("ND1"))
					&& "HIS".equals(resname)) {
				return true;
			} else if ((name.equals("OD1") || name.equals("OD2"))
					&& "ASP".equals(resname)) {
				return true;
			} else if ((name.equals("OE1") || name.equals("OE2"))
					&& "GLU".equals(resname)) {
				return true;
			} else {
				return false;
			}
		}
	}

	public void print() {
		System.out.println("ATOM " + id + " " + name + " " + element + " "
				+ Hnumber + " " + getFx() + " " + getFy() + " " + neig + " "
				+ Hneig + " " + type + " " + flag + " " + wMediated);
	}

}
