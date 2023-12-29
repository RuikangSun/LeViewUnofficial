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
 This class defines a molecule.
 */


import java.io.*;
import java.util.*;

public class Molecule {

	String name;
	LinkedList chains;
	LinkedList structures;
	Properties prop;
	LinkedList water;

	public Molecule() {
		chains = new LinkedList();
		structures = new LinkedList();
		water = new LinkedList();
	}

	public String getName() {
		return name;
	}

	public LinkedList getChains() {
		return chains;
	}

	public Properties getProp() {
		return prop;
	}

	public void setProp(Properties prop) {
		this.prop = prop;
	}

	public LinkedList getWater() {
		return water;
	}

	public void addWater(Atom at) {
		water.add(at);
	}

	public void setName(String c) {
		name = c;
	}

	public void addChain(Chain c) {
		chains.add(c);
	}

	public LinkedList getStruct() {
		return structures;
	}

	public void addStruct(Structure s) {
		structures.add(s);
	}

	public boolean containChain(char c) {
		for (int i = 0; i < chains.size(); i++) {
			Chain chain = (Chain) chains.get(i);
			if (c == chain.getName())
				return true;
		}
		return false;
	}

	public boolean addResidue(int id, Residue r) {
		for (int i = 0; i < chains.size(); i++) {
			Chain chain = (Chain) chains.get(i);
			if (chain.getId() == id) {
				chain.addResidue(r);
				return true;
			}
		}
		return false;
	}

	public boolean addAtom(Chain c, Residue r, Atom a) {
		for (int i = 0; i < chains.size(); i++) {
			Chain chain = (Chain) chains.get(i);
			if (chain.getName() == c.getName()) {
				LinkedList l = chain.getResidues();
				for (int j = 0; j < l.size(); j++) {
					Residue res = (Residue) l.get(j);
					if (res.getId() == r.getId()) {
						res.addAtom(a);
						return true;
					}
				}

			}
		}
		return false;
	}

	public Chain getChain(int id) {
		for (int i = 0; i < chains.size(); i++) {
			Chain chain = (Chain) chains.get(i);
			if (chain.getId() == id)
				return chain;
		}
		return null;
	}

	public void concat(int i1, int i2) {
		Chain c2 = getChain(i2);
		LinkedList l = c2.getResidues();
		for (int i = 0; i < chains.size(); i++) {
			Chain c1 = (Chain) chains.get(i);
			LinkedList ll = c1.getResidues();
			if (c1.getId() == i1) {
				for (int j = 0; j < l.size(); j++) {
					Residue r = (Residue) l.get(j);
					ll.add(r);
				}
				deleteChain(i2);
				break;
			}
		}
	}

	public void deleteChain(int id) {
		for (int i = 0; i < chains.size(); i++) {
			Chain c = (Chain) chains.get(i);
			if (c.getId() == id) {
				chains.remove(i);
				break;
			}
		}
	}

	public void addAtomNeig(Atom id, LinkedList l) {
		for (int v = 0; v < chains.size(); v++) {
			Chain c = (Chain) chains.get(v);
			LinkedList residues = c.getResidues();
			for (int i = 0; i < residues.size(); i++) {
				Residue res = (Residue) residues.get(i);
				LinkedList atoms = res.getAtoms();
				for (int j = 0; j < atoms.size(); j++) {
					Atom at = (Atom) atoms.get(j);
					if (at.getId() == id.getId()) {
						for (int u = 0; u < l.size(); u++) {
							Atom n = (Atom) l.get(u);
							at.addNeig(n);
						}

					}
				}
			}
		}

	}

	public void addAtomHneig(Atom id, LinkedList l) {
		for (int v = 0; v < chains.size(); v++) {
			Chain c = (Chain) chains.get(v);
			LinkedList residues = c.getResidues();
			for (int i = 0; i < residues.size(); i++) {
				Residue res = (Residue) residues.get(i);
				LinkedList atoms = res.getAtoms();
				for (int j = 0; j < atoms.size(); j++) {
					Atom at = (Atom) atoms.get(j);
					if (at.getId() == id.getId()) {
						for (int u = 0; u < l.size(); u++) {
							Atom n = (Atom) l.get(u);
							at.addHneig(n);
						}

					}
				}
			}
		}

	}

	public double getDistance(int id1, int id2) {
		Atom atom1 = new Atom();
		Atom atom2 = new Atom();
		boolean f1 = false;
		boolean f2 = false;
		for (int v = 0; v < chains.size(); v++) {
			Chain c = (Chain) chains.get(v);
			LinkedList residues = c.getResidues();
			for (int i = 0; i < residues.size(); i++) {
				Residue res = (Residue) residues.get(i);
				LinkedList atoms = res.getAtoms();
				for (int j = 0; j < atoms.size(); j++) {
					Atom at = (Atom) atoms.get(j);
					if (at.getId() == id1) {
						atom1 = at;
						f1 = true;
					}
					if (at.getId() == id2) {
						atom2 = at;
						f2 = true;
					}
					if (f1 && f2) {
						return atom1.distance(atom2);
					}
				}
			}
		}
		return -1;
	}

	public double getDistanceRadii(int id1, int id2) {
		Atom atom1 = new Atom();
		Atom atom2 = new Atom();
		boolean f1 = false;
		boolean f2 = false;
		for (int v = 0; v < chains.size(); v++) {
			Chain c = (Chain) chains.get(v);
			LinkedList residues = c.getResidues();
			for (int i = 0; i < residues.size(); i++) {
				Residue res = (Residue) residues.get(i);
				LinkedList atoms = res.getAtoms();
				for (int j = 0; j < atoms.size(); j++) {
					Atom at = (Atom) atoms.get(j);
					if (at.getId() == id1) {
						atom1 = at;
						f1 = true;
					}
					if (at.getId() == id2) {
						atom2 = at;
						f2 = true;
					}
					if (f1 && f2) {

						double result = getRadii(atom1.getElement())
								+ getRadii(atom2.getElement()) + 0.45;
						return result;
					}
				}
			}
		}
		return -1;
	}

	public double getRadii(String e) {

		LinkedList elements = prop.getElements();
		LinkedList radii = prop.getRadii();

		for (int i = 0; i < elements.size(); i++) {
			String el = (String) elements.get(i);
			if (e.equals(el)) {

				double res = (Double) radii.get(i);

				return res;
			}

		}

		return -1;
	}

	public void updateAtomNeig(Atom id) {
		System.out.println("updateAtomNeig");
		for (int v = 0; v < chains.size(); v++) {
			Chain c = (Chain) chains.get(v);
			LinkedList residues = c.getResidues();
			for (int i = 0; i < residues.size(); i++) {
				Residue res = (Residue) residues.get(i);
				LinkedList atoms = res.getAtoms();
				for (int j = 0; j < atoms.size(); j++) {
					Atom at = (Atom) atoms.get(j);
					if (at.getId() == id.getId()) {
						for (int s = 0; s < chains.size(); s++) {
							System.out.println("on y est");
							Chain c2 = (Chain) chains.get(s);
							LinkedList residues2 = c2.getResidues();
							for (int t = 0; t < residues2.size(); t++) {
								Residue res2 = (Residue) residues2.get(t);
								LinkedList atoms2 = res2.getAtoms();
								for (int r = 0; r < atoms2.size(); r++) {
									System.out.println("atom");
									Atom at2 = (Atom) atoms2.get(r);
									if (at2.getId() != id.getId()) {
										double dist = at.distance(at2);
										if ((dist >= 0.4) || (dist <= 1.9)) {
											LinkedList l = new LinkedList();
											l.add(at2);
											addAtomNeig(id, l);
										}
									}

								}
							}
						}

					}
				}
			}
		}

	}

	public void addResNeig(String id1, char c1, String id2, char c2) {

		for (int i = 0; i < chains.size(); i++) {
			Chain chain1 = (Chain) chains.get(i);
			if (chain1.getName() == c1) {
				LinkedList residues1 = chain1.getResidues();
				for (int j = 0; j < residues1.size(); j++) {
					Residue res1 = (Residue) residues1.get(j);
					if (res1.getId().equals(id1)) {

						for (int u = 0; u < chains.size(); u++) {
							Chain chain2 = (Chain) chains.get(u);
							if (chain2.getName() == c2) {
								LinkedList residues2 = chain2.getResidues();
								for (int v = 0; v < residues2.size(); v++) {
									Residue res2 = (Residue) residues2.get(v);
									if (res2.getId().equals(id2)) {
										res1.addResNeig(id2);
									}

								}

							}
						}

					}

				}

			}
		}

	}

	public Atom getAtom(int id) {
		for (int i = 0; i < chains.size(); i++) {
			Chain c = (Chain) chains.get(i);
			LinkedList residues = c.getResidues();
			for (int j = 0; j < residues.size(); j++) {
				Residue res = (Residue) residues.get(j);
				LinkedList atoms = res.getAtoms();
				for (int u = 0; u < atoms.size(); u++) {
					Atom atom = (Atom) atoms.get(u);
					if (atom.getId() == id)
						return atom;
				}
			}

		}
		return null;
	}

	public void print() {
		System.out.println("MOLECULE " + name);
		System.out.println("number of chains " + chains.size());
		System.out.println("number of H2O " + water.size());
		for (int i = 0; i < chains.size(); i++) {
			Chain c = (Chain) chains.get(i);
			c.print();
		}
		System.out
				.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
		for (int i = 0; i < structures.size(); i++) {
			Structure c = (Structure) structures.get(i);
			c.print();
		}

	}

	public void print2() {
		System.out.println("MOLECULE " + name);
		System.out.println("number of chains " + chains.size());
		for (int i = 0; i < chains.size(); i++) {
			Chain c = (Chain) chains.get(i);
			c.print2();
		}

	}

	public void print3() {
		int count = 1;
		for (int i = 0; i < chains.size(); i++) {
			Chain c = (Chain) chains.get(i);
			if (c.getType().equals("protein")) {
				System.out.println("P" + c.getName());
			}
			if (c.getType().equals("dna")) {
				System.out.println("N" + c.getName());
			}

			if (c.getType().equals("ligand")) {
				String line1 = "";
				String line2 = "";
				line1 = line1 + "L00" + count;
				line2 = line2 + "l00" + count + " ";
				count++;
				LinkedList residues = c.getResidues();
				for (int j = 0; j < residues.size(); j++) {
					Residue r = (Residue) residues.get(j);
					if (j == 0) {
						line2 = line2 + r.getName() + " " + r.getId() + "("
								+ c.getName() + ")";
						if (residues.size() != 1) {
							line2 = line2 + " to ";
						}
						line1 = line1 + " " + r.getName() + r.getId() + ":"
								+ c.getName();
					} else {
						if (j == residues.size() - 1) {
							line2 = line2 + r.getName() + " " + r.getId() + "("
									+ c.getName() + ")";
						}
						line1 = line1 + ", " + r.getName() + r.getId() + ":"
								+ c.getName();
					}

				}
				System.out.println(line1);
				System.out.println(line2);
			}
			if (c.getType().equals("ion")) {
				Residue r = (Residue) c.getResidues().getFirst();
				System.out.println("M" + r.getName() + " cyan");
			}
		}

	}

}
