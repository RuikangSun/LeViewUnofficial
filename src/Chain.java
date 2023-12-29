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
 This class defines a chain i.e. protein, nucleic acid, ligand or ion.
 */

import java.io.*;
import java.util.*;

public class Chain {

	char name;
	/* the list of residues contained in this chain */
	LinkedList residues;
	int id;
	boolean path;
	String type;
	Molecule parent;

	public Chain() {
		residues = new LinkedList();
	}

	public Molecule getParent() {
		return parent;
	}

	public void setParent(Molecule m) {
		parent = m;
	}

	public Chain copy() {
		Chain c = new Chain();

		c.setParent(parent);
		c.setId(id);
		c.setName(name);
		c.setType(type);

		LinkedList l = new LinkedList();

		for (int i = 0; i < residues.size(); i++) {
			Residue r = (Residue) residues.get(i);
			l.add(r.copy());
		}
		c.setResidues(l);

		return c;

	}

	public char getName() {
		return name;
	}

	public String getLongName() {
		String tmp = "";
		tmp = tmp + parent.getName() + "_";
		for (int i = 0; i < residues.size(); i++) {
			Residue r = (Residue) residues.get(i);

			tmp = tmp + r.getName() + "_";
		}
		if (name != '?') {
			tmp = tmp + name;
		} else {
			tmp = tmp + id;
		}
		return tmp;
	}

	public String getLongName2() {
		String tmp = "";
		// tmp=tmp+parent.getName()+"_";
		for (int i = 0; i < residues.size(); i++) {
			Residue r = (Residue) residues.get(i);

			tmp = tmp + r.getName() + "_";
		}
		if (name != '?') {
			tmp = tmp + name;
		} else {
			tmp = tmp + id;
		}
		return tmp;
	}

	public LinkedList getResidues() {
		return residues;
	}

	public int getId() {
		return id;
	}

	/* the type can be: protein, nucleic, ligand or ion */
	public String getType() {
		return type;
	}

	public void setName(char c) {
		name = c;
	}

	public void addResidue(Residue r) {
		residues.add(r);
	}

	public void setResidues(LinkedList l) {
		residues = l;
	}

	public void setId(int i) {
		id = i;
	}

	public boolean containRes(String id) {
		for (int i = 0; i < residues.size(); i++) {
			Residue r = (Residue) residues.get(i);
			if (r.getId().equals(id))
				return true;
		}
		return false;
	}

	public int getHetero() {
		int c = 0;
		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);
			c = c + res.getHetero();
		}
		return c;
	}

	/* return the number of atoms contained in this chain */
	public int getNatoms() {
		int c = 0;
		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);
			c = c + res.getNatom();
		}
		return c;
	}

	public void setType(String s) {
		type = s;
	}

	/* return the number of standard and non-standard amino-acids */
	public int getAA() {
		int c = 0;
		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);
			String type = res.getType();
			if ((type.equals("aa")) || (type.equals("ns")))
				c++;
		}
		return c;
	}

	/* return the number of nuclic acids */
	public int getDNA() {
		int c = 0;
		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);
			String type = res.getType();
			if (type.equals("dna"))
				c++;
		}
		return c;
	}

	/* return the number of residues */
	public int getNres() {
		return residues.size();
	}

	/* return true if exist between 2 residues - rescursive algorithm */
	public boolean existPath(Residue res1, Residue res2) {
		path = false;
		for (int i = 0; i < residues.size(); i++) {
			Residue r = (Residue) residues.get(i);
			String u = res1.getId();
			if (r.getId().equals(u)) {
				LinkedList open = new LinkedList();
				LinkedList l = r.getNeig();
				for (int j = 0; j < l.size(); j++) {
					open.add((String) l.get(j));
				}
				LinkedList traite = new LinkedList();
				traite.add(u);
				trait(open, res2.getId(), traite);
				break;
			}
		}
		return path;
	}

	public void trait(LinkedList neig, String end, LinkedList traite) {
		LinkedList temp = new LinkedList();
		if (estDedans(end, neig)) {
			path = true;
		} else {
			if (neig.size() != 0) {
				do {
					if (path == false) {
						String n = (String) neig.get(0);
						neig.remove(0);
						traite.add(n);
						temp.clear();

						for (int i = 0; i < residues.size(); i++) {
							Residue r = (Residue) residues.get(i);

							if (r.getId().equals(n)) {
								LinkedList l = r.getNeig();
								for (int j = 0; j < l.size(); j++) {
									if (estDedans((String) l.get(j), traite)) {
									} else {
										temp.add((String) l.get(j));
									}
								}
								trait(temp, end, traite);
								traite.remove(traite.size() - 1);
								break;
							}
						}
					}
				} while ((path == false) && (neig.size() != 0));
			}
		}
	}

	/* return true if the string e is in the list l */
	public boolean estDedans(String e, LinkedList l) {
		for (int i = 0; i < l.size(); i++) {
			if (e.equals((String) l.get(i)))
				return true;
		}
		return false;
	}

	/* return true is all the residue of the chain are connected to residue r */
	public boolean allConect(Residue r, LinkedList l) {

		for (int i = 0; i < l.size(); i++) {
			Residue res2 = (Residue) l.get(i);
			if (!existPath(r, res2))
				return false;
		}
		return true;
	}

	public void print() {
		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++");
		System.out.println("CHAIN " + id + " " + name);
		System.out.println("number of residues " + residues.size());
		for (int i = 0; i < residues.size(); i++) {
			Residue c = (Residue) residues.get(i);
			c.print();
		}
		System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++");
	}

	public void print2() {
		System.out.println("CHAIN " + id + " " + name + " " + type);
	}

}
