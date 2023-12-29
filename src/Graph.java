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
 This class builds a graph from a ligand and CIF File, makes the mapping, and uptade the ligand
 bond orders.
 */

import java.io.*;
import java.util.*;
import java.net.URL;

public class Graph {

	int[] ids;
	String[] atoms;
	int[][] adja;
	int length;
	LinkedList[] coord;
	String[] elements;
	public LinkedList[][] path;
	LinkedList aaa;
	LinkedList result;

	LinkedList depart;
	LinkedList fin;

	LinkedList ultime;

	Ligand ligand;

	public Graph(Residue res) {
		length = res.getAtoms().size();
		atoms = new String[length];
		adja = new int[length][length];
		ids = new int[length];
		for (int i = 0; i < adja.length; i++) {
			for (int j = 0; j < adja.length; j++) {
				adja[i][j] = 0;
			}
		}
		coord = new LinkedList[length];
		for (int i = 0; i < adja.length; i++) {
			coord[i] = new LinkedList();
		}
		init(res);
		if (atoms != null) {
			traitNames();
		}
	}

	public Graph(String file) {

		init(file);
		if (atoms != null) {
			traitNames();
		}
	}

	public Graph(String[] atoms, int[][] adja, LinkedList[] coord) {
		this.atoms = atoms;
		this.adja = adja;
		this.coord = coord;
		length = atoms.length;

	}

	public Graph(Ligand lig) {
		ligand = lig;
		depart = new LinkedList();
		aaa = new LinkedList();
		fin = new LinkedList();
		result = new LinkedList();
		ultime = new LinkedList();
		init(lig);
		getPath();
		for (int i = 0; i < result.size(); i++) {

			LinkedList l = (LinkedList) result.get(i);
			Atom a1 = (Atom) l.get(0);
			Atom a2 = (Atom) l.get(l.size() - 1);
			LinkedList tmp = new LinkedList();
			for (int j = 1; j < l.size() - 1; j++) {
				Atom at = (Atom) l.get(j);
				tmp.add(at);
			}
			Path p = new Path(a1, a2, tmp);
			ultime.add(p);
		}
	}

	public boolean isUnique() {

		for (int i = 0; i < atoms.length - 1; i++) {

			for (int j = i + 1; j < atoms.length; j++) {
				if (atoms[i].equals(atoms[j]))
					return false;
			}

		}
		return true;
	}

	public LinkedList getUltime() {
		return ultime;
	}

	public void getPath() {

		for (int i = 0; i < depart.size(); i++) {
			for (int j = 0; j < fin.size(); j++) {
				lance((Integer) depart.get(i), (Integer) fin.get(j));
			}

		}

	}

	public void traitNames() {
		for (int i = 0; i < atoms.length; i++) {
			String val = atoms[i];
			if (val.indexOf("\"") != -1) {
				val = val.substring(1, val.length() - 1);
			}
			atoms[i] = val;
		}
	}

	public void lance(int dep, int end) {
		LinkedList open = new LinkedList();
		LinkedList closed = new LinkedList();
		closed.add(dep);

		for (int i = 0; i < adja.length; i++) {
			if (adja[dep][i] == 1)

				open.add(i);
		}

		trait(closed, open, dep, end);
	}

	public void trait(LinkedList closed, LinkedList open, int dep, int end) {

		LinkedList temp = new LinkedList();

		if (open.size() == 0) {
		} else {
			do {
				int element = (Integer) open.get(0);
				Atom a = (Atom) aaa.get(element);
				if (element != end) {
					open.remove(0);
					temp.clear();
					for (int i = 0; i < adja.length; i++) {
						if (adja[element][i] == 1) {

							temp.add(i);
						}

					}
					for (int i = 0; i < temp.size(); i++) {
						if (estDedans((Integer) temp.get(i), closed)) {
							temp.remove(i);
							i--;
						}
					}

					closed.add(element);
					if (closed.size() < 6) {
						trait(closed, temp, dep, end);

					}

				} else {
					open.remove(0);
					closed.add(element);
					LinkedList tt = new LinkedList();
					for (int i = 0; i < closed.size(); i++) {
						int index = (Integer) closed.get(i);
						Atom at = (Atom) aaa.get(index);
						tt.add(at);
					}
					result.add(tt);

				}
				closed.remove(closed.size() - 1);
			} while (open.size() != 0);

		}

	}

	public boolean estDedans(int e, LinkedList l) {
		for (int i = 0; i < l.size(); i++) {
			if (e == (Integer) l.get(i))
				return true;
		}
		return false;
	}

	public void init(Residue res) {
		LinkedList l = res.getAtoms();

		for (int i = 0; i < l.size(); i++) {
			Atom at = (Atom) l.get(i);
			ids[i] = at.getId();

			atoms[i] = at.getName();
			coord[i].add(at.getX());
			coord[i].add(at.getY());
			coord[i].add(at.getZ());
			LinkedList temp = at.getNeig();

			for (int j = 0; j < temp.size(); j++) {
				Atom att = (Atom) temp.get(j);
				int ID = att.getId();
				int t = res.getAtomPos(ID);

				if (t != -1) {
					adja[i][t] = 1;
				}
			}
		}

	}

	public boolean isInside(Atom a, LinkedList l) {

		for (int i = 0; i < l.size(); i++) {
			Atom at = (Atom) l.get(i);
			if (at.getId() == a.getId())
				return true;
		}
		return false;
	}

	public void init(Ligand lig) {

		LinkedList w = lig.getWater();

		LinkedList l = new LinkedList();

		for (int i = 0; i < w.size(); i++) {
			Atom at = (Atom) w.get(i);

			LinkedList l1 = at.getwMediated1();
			LinkedList l2 = at.getwMediated2();
			LinkedList l3 = at.getwMediated3();

			if (l1.size() != 0) {
				if (!isInside(at, l)) {
					l.add(at);
				}
				for (int j = 0; j < l1.size(); j++) {
					Atom at2 = (Atom) l1.get(j);
					at2.setwType(1);
					if (!isInside(at2, l)) {
						l.add(at2);
					}

				}
			}

			if (l2.size() != 0) {
				if (!isInside(at, l)) {
					l.add(at);
				}
				for (int j = 0; j < l2.size(); j++) {
					Atom at2 = (Atom) l2.get(j);
					at2.setwType(2);
					if (!isInside(at2, l)) {
						l.add(at2);

					}
				}
			}

			if (l3.size() != 0) {
				if (!isInside(at, l)) {
					l.add(at);
				}
				for (int j = 0; j < l3.size(); j++) {
					Atom at2 = (Atom) l3.get(j);
					at2.setwType(3);
					if (!isInside(at2, l)) {
						l.add(at2);
					}

				}
			}

		}

		LinkedList depAt = new LinkedList();
		LinkedList finAt = new LinkedList();

		for (int i = 0; i < l.size(); i++) {
			Atom at = (Atom) l.get(i);

			if (at.getwType() == 1) {
				depart.add(i);
				depAt.add(at);
			}
			if (at.getwType() == 3) {
				fin.add(i);
				finAt.add(at);
			}

		}

		length = l.size();
		ids = new int[length];
		atoms = new String[length];
		adja = new int[length][length];
		aaa = l;

		for (int i = 0; i < l.size(); i++) {
			Atom at = (Atom) l.get(i);
			ids[i] = at.getId();
			atoms[i] = at.getName();
		}

		for (int i = 0; i < adja.length; i++) {
			for (int j = 0; j < adja.length; j++) {
				adja[i][j] = 0;
			}
		}

		/* adjacency matrix update */
		for (int i = 0; i < l.size(); i++) {

			Atom at = (Atom) l.get(i);

			LinkedList l1 = at.getwMediated1();

			for (int j = 0; j < l1.size(); j++) {
				Atom at2 = (Atom) l1.get(j);
				adja[getIndex(at)][getIndex(at2)] = 1;
			}

			LinkedList l2 = at.getwMediated2();

			for (int j = 0; j < l2.size(); j++) {
				Atom at2 = (Atom) l2.get(j);
				adja[getIndex(at)][getIndex(at2)] = 1;
			}

			LinkedList l3 = at.getwMediated3();

			for (int j = 0; j < l3.size(); j++) {
				Atom at2 = (Atom) l3.get(j);
				adja[getIndex(at)][getIndex(at2)] = 1;
			}

		}

	}

	public int getIndex(Atom at) {
		int id = at.getId();

		for (int i = 0; i < ids.length; i++) {
			if (ids[i] == id)
				return i;
		}
		return -1;
	}

	public int getIndex(Atom at, LinkedList l) {
		int id = at.getId();

		for (int i = 0; i < l.size(); i++) {
			Atom a2 = (Atom) l.get(i);
			if (a2.getId() == id)
				return i;
		}
		return -1;
	}

	public void init(String file) {
		String fic = "";
		fic = "mmcif/" + file + ".cif";

		LinkedList temp = new LinkedList();
		LinkedList temp2 = new LinkedList();
		int atomCount = 0;
		int bondCount = 0;
		try {

			URL monUrl = getClass().getResource(fic);

			InputStreamReader ins = new InputStreamReader(monUrl.openStream());

			BufferedReader buff = new BufferedReader(ins);

			try {

				String line;

				boolean fin = false;
				boolean deb = false;
				boolean deb2 = false;
				boolean created = false;

				boolean finB = false;
				boolean debB = false;
				boolean deb2B = false;

				while ((line = buff.readLine()) != null) {

					/* atoms */

					if(line.length()!=0){
					//System.out.println(line.length());
					if (line.indexOf("chem_comp_atom") != -1) {
						deb = true;
						line = buff.readLine();
					}
					if ((deb == true)) {
						if (line.charAt(0) != '_') {
							deb2 = true;
						}
						if (line.indexOf("#") != -1) {
							fin = true;
						}
					}

					if ((deb == true) && (deb2 == true) && (fin == false)) {

						StringTokenizer tok = new StringTokenizer(line, " ");
						tok.nextToken();
						String id = tok.nextToken();
						tok.nextToken();
						String e = tok.nextToken();
						temp.add(id);
						temp2.add(e);

					}

					if ((deb == true) && (deb2 == true) && (fin == true)
							&& (created == false)) {
						atoms = new String[temp.size()];
						elements = new String[temp.size()];
						length = temp.size();
						atomCount = atoms.length;
						for (int i = 0; i < temp.size(); i++) {

							atoms[i] = (String) temp.get(i);
							elements[i] = (String) temp2.get(i);

						}

						adja = new int[length][length];

						for (int i = 0; i < adja.length; i++) {
							for (int j = 0; j < adja.length; j++) {
								adja[i][j] = 0;
							}
						}

						created = true;
					}

					/* Bonds */

					if (line.indexOf("chem_comp_bond.comp") != -1) {
						debB = true;
						line = buff.readLine();
					}
					if ((debB == true)) {
						if (line.charAt(0) != '_') {
							deb2B = true;
						}
						if (line.indexOf("#") != -1) {
							finB = true;
						}
					}

					if ((debB == true) && (deb2B == true) && (finB == false)) {

						StringTokenizer tok = new StringTokenizer(line, " ");
						tok.nextToken();
						String id1 = tok.nextToken();
						String id2 = tok.nextToken();
						String order = tok.nextToken();

						adja[getPos(id1)][getPos(id2)] = getOrder(order);
						adja[getPos(id2)][getPos(id1)] = getOrder(order);

					}
					}

				}

			} finally {
				buff.close();
			}
		} catch (IOException ioe) {
			System.out.println("Error --" + ioe.toString());
		}

	}

	public boolean isCompa(String[] t) {

		for (int i = 0; i < t.length; i++) {
			if (getPos(t[i]) == -1)
				return false;
		}
		return true;
	}

	public int getPos(String s) {
		for (int i = 0; i < atoms.length; i++) {
			if (atoms[i].equals(s))
				return i;
		}
		return -1;
	}

	public int getOrder(String s) {
		if (s.equals("SING"))
			return 1;
		if (s.equals("DOUB"))
			return 2;
		if (s.equals("TRIP"))
			return 3;
		if (s.equals("QUAD"))
			return 4;
		if (s.equals("AROM"))
			return 6;
		return -1;
	}

	public void print() {
		for (int i = 0; i < length; i++) {
			System.out.print(atoms[i] + " ; ");
		}
		System.out.println("");
		for (int i = 0; i < adja.length; i++) {
			for (int j = 0; j < adja.length; j++) {
				System.out.print(adja[i][j] + " | ");
			}
			System.out.println("");
		}
	}

	public Graph withoutH() {
		LinkedList l = new LinkedList();
		LinkedList pos = new LinkedList();
		for (int i = 0; i < length; i++) {
			if (!elements[i].equals("H")) {
				l.add(atoms[i]);
				pos.add(i);
			}
		}

		String[] at = new String[l.size()];
		int[][] adj = new int[l.size()][l.size()];
		LinkedList[] coo = new LinkedList[l.size()];

		for (int i = 0; i < l.size(); i++) {
			String temp = (String) l.get(i);
			int pp = (Integer) pos.get(i);
			at[i] = temp;

			for (int j = 0; j < l.size(); j++) {
				adj[i][j] = adja[(Integer) pos.get(i)][(Integer) pos.get(j)];
			}

		}

		Graph g = new Graph(at, adj, coo);
		return g;

	}

	public int getH(int pos) {
		int res = 0;
		for (int i = 0; i < length; i++) {
			if (adja[pos][i] != 0) {
				if (elements[i].equals("H"))
					res++;
			}
		}
		return res;
	}

}
