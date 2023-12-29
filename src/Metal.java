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
 This class defines a metal.
 */

import java.io.*;
import java.util.*;

public class Metal {

	Atom atom;
	Chain chain;
	LinkedList close;
	Residue residue;
	Molecule mol;
	/* contains explicit residues */
	LinkedList resExplicit;
	LinkedList explicit;
	LinkedList explicit2;
	/* contains ligand chains */
	LinkedList ligands;
	LinkedList atoms;
	LinkedList lightLigands;
	LinkedList waterAll;
	LinkedList water;
	LinkedList contactList;

	LinkedList posName;
	LinkedList posName2;

	double cutoff;

	public Metal() {

	}

	public Metal(Chain c, double c1) {

		cutoff = c1;
		chain = c;
		mol = chain.getParent();
		LinkedList res = chain.getResidues();
		residue = (Residue) res.get(0);
		LinkedList a = residue.getAtoms();
		lightLigands = new LinkedList();
		atom = (Atom) a.get(0);

		atom.setFx(atom.getX());
		atom.setFy(atom.getY());

		atoms = new LinkedList();
		atoms.add(atom);

		explicit = new LinkedList();
		explicit2 = new LinkedList();

		updateExplicit();

		resExplicit = new LinkedList();
		ligands = new LinkedList();
		waterAll = mol.getWater();
		water = new LinkedList();
		contactList = new LinkedList();

		posName = new LinkedList();
		posName2 = new LinkedList();

		close = getClose();

		updateLigands();

		traitLigands();

		getResExplicit();
		getWater();

		updateClose();

		contactList = getContactList();
		traitContactList();

		posName = getPoseName();
		posName2 = getPoseName2();

	}

	public void updateLigands() {
		LinkedList n = atom.getNeig();
		for (int i = 0; i < n.size(); i++) {
			Atom at = (Atom) n.get(i);
			Residue rr = at.getParent();
			Chain c = rr.getParent();
			if (c.getType().equals("ligand")) {
				if (!isInsideLig(rr)) {
					ligands.add(rr.getParent());
				}
			}
		}

		LinkedList hn = atom.getHneig();
		for (int i = 0; i < hn.size(); i++) {
			Atom at = (Atom) hn.get(i);
			Residue rr = at.getParent();
			Chain c = rr.getParent();
			if (c.getType().equals("ligand")) {
				if (!isInsideLig(rr)) {
					ligands.add(rr.getParent());
				}
			}
		}

	}

	public Chain getParent() {
		return chain;
	}

	public Molecule getMol() {
		return mol;
	}

	public LinkedList getLightLigands() {
		return lightLigands;
	}

	public LinkedList getExplicit() {
		return explicit;
	}

	public LinkedList getExplicit2() {
		return explicit2;
	}

	public Atom getAtom() {
		return atom;
	}

	public LinkedList getPoseName() {

		LinkedList res = new LinkedList();

		for (int i = 0; i < lightLigands.size(); i++) {
			LigandLight lig = (LigandLight) lightLigands.get(i);

			LinkedList rr = lig.getResidues();
			for (int j = 0; j < rr.size(); j++) {
				Residue r = (Residue) rr.get(j);

				res.add(getPosName(r));
			}

		}

		return res;

	}

	public void updateClose() {

		LinkedList copy = new LinkedList();
		for (int i = 0; i < close.size(); i++) {
			Residue r = (Residue) close.get(i);

			if (!isInsideEx(r)) {
				if (!isInsideLig(r)) {
					copy.add(r);
				}

			}

		}

		close = copy;
	}

	public boolean isInsideEx(Residue r) {

		for (int i = 0; i < resExplicit.size(); i++) {
			Explicit e = (Explicit) resExplicit.get(i);
			Residue res = e.getResidue();
			if (res == r)
				return true;

		}
		return false;
	}

	public boolean isInsideLig(Residue r) {

		for (int i = 0; i < ligands.size(); i++) {
			Chain c = (Chain) ligands.get(i);

			if (r.getParent() == c)
				return true;

		}
		return false;
	}

	public LinkedList getPoseName2() {

		LinkedList res = new LinkedList();

		for (int i = 0; i < resExplicit.size(); i++) {
			Explicit e = (Explicit) resExplicit.get(i);
			LinkedList tmp = e.getResidues();
			Residue r = (Residue) tmp.get(0);

			res.add(getPosName2(r));

		}

		return res;

	}

	public boolean conflit(Atom atom, int pos) {

		double r = 2;

		for (int i = 0; i < atoms.size(); i++) {
			Atom a2 = (Atom) atoms.get(i);
			double dd = Math.sqrt((atom.getFx() - a2.getFx())
					* (atom.getFx() - a2.getFx()) + (atom.getFy() - a2.getFy())
					* (atom.getFy() - a2.getFy()));
			if (dd <= r) {

				return true;
			}

		}

		if (resExplicit.size() > 0) {

			for (int u = 0; u < resExplicit.size(); u++) {

				Explicit ee = (Explicit) resExplicit.get(u);

				LinkedList tmp = ee.getAtoms();

				for (int i = 0; i < tmp.size(); i++) {
					Atom a2 = (Atom) tmp.get(i);
					double dd = Math.sqrt((atom.getFx() - a2.getFx())
							* (atom.getFx() - a2.getFx())
							+ (atom.getFy() - a2.getFy())
							* (atom.getFy() - a2.getFy()));
					if (dd <= r) {

						return true;
					}
				}
			}
		}

		for (int i = 0; i < pos; i++) {
			Contact con = (Contact) contactList.get(i);
			Atom a2 = con.getTrait();
			double dd = Math.sqrt((atom.getFx() - a2.getFx())
					* (atom.getFx() - a2.getFx()) + (atom.getFy() - a2.getFy())
					* (atom.getFy() - a2.getFy()));

			if (dd <= r) {

				return true;
			}
		}

		return false;

	}

	public LinkedList getContacts() {
		return contactList;
	}

	public void rotation(Contact con, double angle) {

		angle = angle * Math.PI / 180;

		Atom at = con.getTrait();
		Atom pivot = con.getRef();

		double OAx = at.getFx() - pivot.getFx();
		double OAy = at.getFy() - pivot.getFy();

		double resx = pivot.getFx() + Math.cos(angle) * OAx - Math.sin(angle)
				* OAy;
		double resy = pivot.getFy() + Math.sin(angle) * OAx + Math.cos(angle)
				* OAy;

		con.setX(resx);
		con.setY(resy);

	}

	public double getNumberAtom(Hbond h, LinkedList l, double r) {

		double res = 0;
		int nbre = 0;

		Atom atom = h.getAtomRes();

		for (int i = 0; i < atoms.size(); i++) {
			Atom a2 = (Atom) atoms.get(i);
			double dd = Math.sqrt((atom.getFx() - a2.getFx())
					* (atom.getFx() - a2.getFx()) + (atom.getFy() - a2.getFy())
					* (atom.getFy() - a2.getFy()));

			if (dd <= r) {

				res = res + dd;
				nbre++;
			}
		}

		if (resExplicit.size() > 0) {

			for (int u = 0; u < resExplicit.size(); u++) {

				Explicit ee = (Explicit) resExplicit.get(u);

				LinkedList tmp = ee.getAtoms();

				for (int i = 0; i < tmp.size(); i++) {
					Atom a2 = (Atom) tmp.get(i);
					double dd = Math.sqrt((atom.getFx() - a2.getFx())
							* (atom.getFx() - a2.getFx())
							+ (atom.getFy() - a2.getFy())
							* (atom.getFy() - a2.getFy()));
					if (dd <= r) {

						res = res + dd;
						nbre++;
					}
				}
			}

		}

		for (int i = 0; i < l.size(); i++) {
			Hbond hb = (Hbond) l.get(i);
			Atom a2 = hb.getAtomRes();

			if (h.equals(hb))
				break;

			double dd = Math.sqrt((atom.getFx() - a2.getFx())
					* (atom.getFx() - a2.getFx()) + (atom.getFy() - a2.getFy())
					* (atom.getFy() - a2.getFy()));

			if (dd <= r) {

				res = res + dd;
				nbre++;
			}

		}

		if (nbre != 0) {
			res = res / nbre;
		} else {
			res = 1000;
		}

		return res;
	}

	public LinkedList returnPoseName() {
		return posName;
	}

	public LinkedList returnPoseName2() {
		return posName2;
	}

	public LinkedList getContactList() {
		LinkedList res = new LinkedList();

		for (int i = 0; i < close.size(); i++) {
			Residue r = (Residue) close.get(i);

			Atom ref = null;
			Atom trait = null;

			LinkedList at = r.getAtoms();
			LinkedList l = barycenter2(at);

			double x0 = (Double) l.get(0);
			double y0 = (Double) l.get(1);

			double t = 100000;

			for (int u = 0; u < at.size(); u++) {
				Atom a = (Atom) at.get(u);

				double dd = Math.sqrt((x0 - a.getX()) * (x0 - a.getX())
						+ (y0 - a.getY()) * (y0 - a.getY()));
				if (dd < t) {
					t = dd;
					trait = a;
				}

			}

			double d = 100000;

			for (int u = 0; u < at.size(); u++) {
				Atom a = (Atom) at.get(u);
				for (int j = 0; j < atoms.size(); j++) {
					Atom atom = (Atom) atoms.get(j);
					double dist = a.distance(atom);
					if (dist < d) {
						d = dist;
						ref = atom;
						trait = a;
					}
				}
			}
			Contact con = new Contact(r, trait, ref, d);
			res.add(con);

		}

		return res;

	}

	public void traitContactList() {

		for (int i = 0; i < contactList.size(); i++) {
			Contact con = (Contact) contactList.get(i);

			Atom ref = con.getRef();

			Atom res = ref.clone();
			res.setFx(res.getFx() + 1);
			res.setFy(res.getFy() + 1);

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = res.getFx();
			double y1 = res.getFy();
			double x2 = ref.getFx();
			double y2 = ref.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x2;
			double yo = y2;
			double r = con.getDistance();

			if (b != 0.0) {
				double u = a * (c + b * yo) - b * b * xo;
				double v = a * a + b * b;
				double w = c + b * yo;
				double deltap = u * u - v * (b * b * (xo * xo - r * r) + w * w);
				if (deltap >= 0.0) {
					sol1x = (-u + 1.0 * Math.sqrt(deltap)) / v;
					sol1y = -(a * sol1x + c) / b;
					sol2x = (-u + -1.0 * Math.sqrt(deltap)) / v;
					sol2y = -(a * sol2x + c) / b;
				}

			} else {
				sol1x = -c / a;
				double u = sol1x - xo;
				double v = r * r - u * u;
				if (v >= 0.0) {
					sol1y = yo + 1.0 * Math.sqrt(v);
					sol2y = yo + -1.0 * Math.sqrt(v);
					sol2x = sol1x;
				}

			}

			double d1 = Math.sqrt((ref.getFx() - sol1x) * (ref.getFx() - sol1x)
					+ (ref.getFy() - sol1y) * (ref.getFy() - sol1y));
			double d2 = Math.sqrt((ref.getFx() - sol2x) * (ref.getFx() - sol2x)
					+ (ref.getFy() - sol2y) * (ref.getFy() - sol2y));

			double solx = 0;
			double soly = 0;

			if (d1 < d2) {

				solx = sol2x;
				soly = sol2y;

			} else {

				solx = sol1x;
				soly = sol1y;
			}

			con.setX(solx);
			con.setY(soly);

			if (conflit(con.getTrait(), i)) {

				double[] rot = new double[72];

				for (int s = 0; s < rot.length; s++) {
					rot[s] = 5;
				}

				for (int v = 0; v < rot.length; v++) {

					if (conflit(con.getTrait(), i)) {

						rotation(con, rot[v]);
					} else {
						break;
					}

				}

			}

		}

	}

	public double[] getPosName(Residue r) {

		double[] tab = new double[2];

		LinkedList at = r.getAtoms();

		LinkedList temp = barycenter(at);

		tab[0] = (Double) temp.get(0);
		tab[1] = (Double) temp.get(1);

		double xo = tab[0];
		double yo = tab[1];

		double distMax = r.getDistMax(xo, yo) + 0.5;

		if (getNumberAtom(tab, 2) == 100) {
			return tab;
		}

		double angle = 10;
		int o = 0;
		double tmp = 0;
		double dep = 0;

		tab[0] = tab[0] + 0.5;

		for (int j = 0; j * 0.5 < distMax; j++) {

			for (int i = 0; i < 36; i++) {

				rotation(tab, xo, yo, 10);
				double n = getNumberAtom(tab, 2);

				if (n == 100) {
					return tab;
				} else {
					if (n > tmp) {
						tmp = n;
						o = i;
					}
				}

			}

			tab[0] = tab[0] + 0.5;
			dep = j;

		}

		tab[0] = xo + 0.5 + dep * 0.5;
		tab[1] = yo;

		rotation(tab, xo, yo, 10 + o * 10);

		return tab;

	}

	public double[] getPosName2(Residue r) {

		if (r.getType().equals("ion")) {
			return null;
		}

		double[] tab = new double[2];

		LinkedList at = r.getAtoms();

		LinkedList temp = barycenter(at);

		tab[0] = (Double) temp.get(0);
		tab[1] = (Double) temp.get(1);

		double xo = tab[0];
		double yo = tab[1];

		double distMax = r.getDistMax(xo, yo) + 0.5;

		if (getNumberAtom2(tab, 2) == 100) {
			return tab;
		}

		double angle = 10;
		int o = 0;
		double tmp = 0;
		double dep = 0;

		tab[0] = tab[0] + 0.5;

		for (int j = 0; j * 0.5 < distMax; j++) {

			for (int i = 0; i < 36; i++) {

				rotation(tab, xo, yo, 10);
				double n = getNumberAtom2(tab, 2);

				if (n == 100) {
					return tab;
				} else {
					if (n > tmp) {
						tmp = n;
						o = i;
					}
				}

			}

			tab[0] = tab[0] + 0.5;
			dep = j;

		}

		tab[0] = xo + 0.5 + dep * 0.5;
		tab[1] = yo;

		rotation(tab, xo, yo, 10 + o * 10);

		return tab;

	}

	public double getNumberAtom2(double[] tab, double r) {

		double res = 0;
		int nbre = 0;

		double x = tab[0];
		double y = tab[1];

		for (int i = 0; i < atoms.size(); i++) {
			Atom a2 = (Atom) atoms.get(i);
			double dd = Math.sqrt((x - a2.getFx()) * (x - a2.getFx())
					+ (y - a2.getFy()) * (y - a2.getFy()));

			if (dd <= r) {

				res = res + dd;
				nbre++;
			}
		}

		if (resExplicit.size() > 0) {

			for (int u = 0; u < resExplicit.size(); u++) {

				Explicit ee = (Explicit) resExplicit.get(u);

				LinkedList tmp = ee.getAtoms();

				for (int i = 0; i < tmp.size(); i++) {
					Atom a2 = (Atom) tmp.get(i);
					double dd = Math.sqrt((x - a2.getFx()) * (x - a2.getFx())
							+ (y - a2.getFy()) * (y - a2.getFy()));
					if (dd <= r) {

						res = res + dd;
						nbre++;
					}
				}
			}

		}

		for (int u = 0; u < lightLigands.size(); u++) {
			LigandLight lig = (LigandLight) lightLigands.get(u);
			LinkedList att = lig.getAtoms();

			for (int i = 0; i < att.size(); i++) {
				Atom a2 = (Atom) att.get(i);

				double dd = Math.sqrt((x - a2.getFx()) * (x - a2.getFx())
						+ (y - a2.getFy()) * (y - a2.getFy()));
				if (dd < r) {

					res = res + dd;
					nbre++;
				}
			}

		}

		for (int i = 0; i < water.size(); i++) {
			Atom a2 = (Atom) water.get(i);

			double dd = Math.sqrt((x - a2.getFx()) * (x - a2.getFx())
					+ (y - a2.getFy()) * (y - a2.getFy()));
			if (dd < r) {

				res = res + dd;
				nbre++;
			}
		}

		for (int i = 0; i < posName.size(); i++) {
			double[] tabe = (double[]) posName.get(i);
			double xx = tabe[0];
			double yy = tabe[1];
			double dd = Math.sqrt((x - xx) * (x - xx) + (y - yy) * (y - yy));

			if (dd <= r) {

				res = res + dd;
				nbre++;
			}
		}

		if (nbre != 0) {
			res = res / nbre;
		} else {
			res = 100;
		}

		return res;
	}

	public double getNumberAtom(double[] tab, double r) {

		double res = 0;
		int nbre = 0;

		double x = tab[0];
		double y = tab[1];

		for (int i = 0; i < atoms.size(); i++) {
			Atom a2 = (Atom) atoms.get(i);
			double dd = Math.sqrt((x - a2.getFx()) * (x - a2.getFx())
					+ (y - a2.getFy()) * (y - a2.getFy()));

			if (dd <= r) {

				res = res + dd;
				nbre++;
			}
		}

		if (resExplicit.size() > 0) {

			for (int u = 0; u < resExplicit.size(); u++) {

				Explicit ee = (Explicit) resExplicit.get(u);

				LinkedList tmp = ee.getAtoms();

				for (int i = 0; i < tmp.size(); i++) {
					Atom a2 = (Atom) tmp.get(i);
					double dd = Math.sqrt((x - a2.getFx()) * (x - a2.getFx())
							+ (y - a2.getFy()) * (y - a2.getFy()));
					if (dd <= r) {

						res = res + dd;
						nbre++;
					}
				}
			}

		}

		for (int i = 0; i < water.size(); i++) {
			Atom a2 = (Atom) water.get(i);

			double dd = Math.sqrt((x - a2.getFx()) * (x - a2.getFx())
					+ (y - a2.getFy()) * (y - a2.getFy()));
			if (dd < r) {

				res = res + dd;
				nbre++;
			}
		}

		for (int i = 0; i < posName.size(); i++) {
			double[] tabe = (double[]) posName.get(i);
			double xx = tabe[0];
			double yy = tabe[1];
			double dd = Math.sqrt((x - xx) * (x - xx) + (y - yy) * (y - yy));

			if (dd <= r) {

				res = res + dd;
				nbre++;
			}
		}

		if (nbre != 0) {
			res = res / nbre;
		} else {
			res = 100;
		}

		return res;
	}

	public void rotation(double[] tab, double xo, double yo, double angle) {

		angle = angle * Math.PI / 180;

		double OAx = tab[0] - xo;
		double OAy = tab[1] - yo;

		double resx = xo + Math.cos(angle) * OAx - Math.sin(angle) * OAy;
		double resy = yo + Math.sin(angle) * OAx + Math.cos(angle) * OAy;

		tab[0] = resx;
		tab[1] = resy;

	}

	public LinkedList getClose() {
		LinkedList l = new LinkedList();
		boolean found = false;

		LinkedList chains = mol.getChains();

		for (int c = 0; c < chains.size(); c++) {
			Chain ch = (Chain) chains.get(c);
			if ((ch.getId() != chain.getId()) && (!ch.getType().equals("ion"))) {

				LinkedList residues = ch.getResidues();

				for (int j = 0; j < residues.size(); j++) {
					Residue r = (Residue) residues.get(j);

					if (residue.getId() != r.getId()) {
						found = false;

						LinkedList atoms2 = r.getAtoms();

						if (found)
							break;
						for (int v = 0; v < atoms2.size(); v++) {
							Atom a2 = (Atom) atoms2.get(v);
							if (isDist(atom, a2, 5.0)) {
								if (r.getParent().getType().equals("ligand")) {
									if (!isInside(r.getParent(), l)) {

									}

								} else {

									l.add(r);
									found = true;
									break;
								}

							}

						}

					}

				}
			}
		}
		return l;
	}

	public boolean isDist(Atom a1, Atom a2, double d) {
		if (a1.distance(a2) < d)
			return true;
		else
			return false;

	}

	public void updateExplicit() {

		LinkedList hn = atom.getHneig();
		for (int i = 0; i < hn.size(); i++) {
			Atom at = (Atom) hn.get(i);
			at.init();
			Bond b = new Bond(atom, at, 1);

			explicit.add(b);
		}

		LinkedList n = atom.getNeig();
		for (int i = 0; i < n.size(); i++) {
			Atom at = (Atom) n.get(i);
			at.init();
			Residue r = at.getParent();
			if (!isWater(r.getName())) {
				Bond b = new Bond(atom, at, 1);

				explicit.add(b);
			}
		}

	}

	public boolean isWater(String s) {
		if ((s.equals("HOH")) || (s.equals("H2O")) || (s.equals("WAT")))
			return true;
		else
			return false;
	}

	public boolean isInside(Chain c, LinkedList l) {

		for (int i = 0; i < l.size(); i++) {

			Residue rr = (Residue) l.get(i);
			Chain cc = rr.getParent();

			if (c.getId() == cc.getId())
				return true;

		}
		return false;

	}

	public void traitLigands() {

		for (int i = 0; i < ligands.size(); i++) {
			Chain ch = (Chain) ligands.get(i);
			LigandLight lig = new LigandLight(ch);

			atom.setFx(atom.getX());
			atom.setFy(atom.getY());

			/* ligand atoms linked to metal */
			LinkedList at = lig.getAtoms();
			LinkedList lie = new LinkedList();
			for (int u = 0; u < at.size(); u++) {
				Atom a = (Atom) at.get(u);
				LinkedList hn = a.getHneig();
				for (int v = 0; v < hn.size(); v++) {
					Atom a2 = (Atom) hn.get(v);
					if (a2.getId() == atom.getId())
						lie.add(a);

				}
			}

			Atom res = new Atom();

			if (lie.size() == 1) {
				res = (Atom) lie.get(0);
			} else {

				LinkedList tmp = barycenter(lie);
				res.setFx((Double) tmp.get(0));
				res.setFy((Double) tmp.get(1));
			}

			double dist = 0;

			for (int u = 0; u < lie.size(); u++) {
				Atom aaa = (Atom) lie.get(u);

				dist = dist + atom.distance(aaa);

			}

			dist = dist / lie.size();

			Atom ref = atom;

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = res.getFx();
			double y1 = res.getFy();
			double x2 = ref.getFx();
			double y2 = ref.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = dist;

			if (b != 0.0) {
				double u = a * (c + b * yo) - b * b * xo;
				double v = a * a + b * b;
				double w = c + b * yo;
				double deltap = u * u - v * (b * b * (xo * xo - r * r) + w * w);
				if (deltap >= 0.0) {
					sol1x = (-u + 1.0 * Math.sqrt(deltap)) / v;
					sol1y = -(a * sol1x + c) / b;
					sol2x = (-u + -1.0 * Math.sqrt(deltap)) / v;
					sol2y = -(a * sol2x + c) / b;
				}

			} else {
				sol1x = -c / a;
				double u = sol1x - xo;
				double v = r * r - u * u;
				if (v >= 0.0) {
					sol1y = yo + 1.0 * Math.sqrt(v);
					sol2y = yo + -1.0 * Math.sqrt(v);
					sol2x = sol1x;
				}

			}

			double d1 = Math.sqrt((ref.getFx() - sol1x) * (ref.getFx() - sol1x)
					+ (ref.getFy() - sol1y) * (ref.getFy() - sol1y));
			double d2 = Math.sqrt((ref.getFx() - sol2x) * (ref.getFx() - sol2x)
					+ (ref.getFy() - sol2y) * (ref.getFy() - sol2y));

			double solx = 0;
			double soly = 0;

			if (d1 < d2) {

				solx = sol2x;
				soly = sol2y;

			} else {

				solx = sol1x;
				soly = sol1y;
			}

			double ddx, ddy;

			ddx = solx - res.getFx();
			ddy = soly - res.getFy();

			lig.translation(ddx, ddy);

			double[] rot = new double[36];

			for (int s = 0; s < rot.length; s++) {
				rot[s] = 10;
			}

			boolean bbb = false;

			double dist2 = 5.0;
			double nb = getNumberInter(lig);

			double na = getNumberAtom(lig, dist2);
			double op1 = 0;
			double op2 = 0;

			for (int v = 0; v < rot.length; v++) {

				lig.rotation(atom, rot[v]);

				for (int u = 0; u < rot.length; u++) {

					lig.rotation(res, rot[u]);

					double nb2 = getNumberInter(lig);

					double na2 = getNumberAtom(lig, dist2);

					if (nb2 < nb) {
						nb = nb2;
						na = na2;
						op1 = v;
						op2 = u;
					} else {
						if (nb2 == nb) {
							if (na2 > na) {
								nb = nb2;
								na = na2;
								op1 = v;
								op2 = u;
							}

						}

					}

				}

			}

			lig.rotation(atom, (op1 + 1) * 10);
			lig.rotation(res, (op2 + 1) * 10);

			lightLigands.add(lig);
		}

	}

	public void getResExplicit() {

		LinkedList res = new LinkedList();
		LinkedList aaa = new LinkedList();

		for (int i = 0; i < explicit.size(); i++) {
			Bond b = (Bond) explicit.get(i);

			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();

			Residue r1 = a1.getParent();
			Residue r2 = a2.getParent();

			if ((!r1.getParent().getType().equals("ligand"))
					&& (!r2.getParent().getType().equals("ligand"))) {
				if (a1.getId() != atom.getId()) {
					res.add(r1);
					aaa.add(a1);
				} else {
					res.add(r2);
					aaa.add(a2);
				}

			}
		}

		for (int i = 0; i < res.size(); i++) {
			Atom aa = (Atom) aaa.get(i);
			Residue rrr = (Residue) res.get(i);

			Explicit e = new Explicit(rrr, aa);

			traitExplicit(e, aa);

			double[] rot = new double[36];

			for (int s = 0; s < rot.length; s++) {
				rot[s] = 10;
			}

			double dist2 = 5.0;
			double nb = getNumberInter(e);
			double na = getNumberAtom(e, dist2);
			double op1 = 0;
			double op2 = 0;

			for (int v = 0; v < rot.length; v++) {

				e.rotation(aa, rot[v]);

				for (int u = 0; u < rot.length; u++) {

					e.rotation(e.getAtLigand(), rot[u]);

					double nb2 = getNumberInter(e);
					double na2 = getNumberAtom(e, dist2);

					if (nb2 < nb) {
						nb = nb2;
						na = na2;
						op1 = v;
						op2 = u;
					} else {
						if (nb2 == nb) {
							if (na2 > na) {
								nb = nb2;
								na = na2;
								op1 = v;
								op2 = u;
							}

						}

					}

				}

			}

			e.rotation(aa, (op1 + 1) * 10);
			e.rotation(e.getAtLigand(), (op2 + 1) * 10);

			resExplicit.add(e);
		}

	}

	public LinkedList getWaterList() {
		return water;
	}

	public void getWater() {

		for (int i = 0; i < waterAll.size(); i++) {
			Atom w = (Atom) waterAll.get(i);
			double dist = w.distance(atom);
			if (dist < 4.0) {

				w.setFx(w.getX());
				w.setFy(w.getY());

				double[] rot = new double[36];

				for (int s = 0; s < rot.length; s++) {
					rot[s] = 10;
				}

				double dist2 = 5.0;
				double nb = getNumberInter(w);
				double na = getNumberAtom(w, dist2);
				double op1 = 0;

				for (int v = 0; v < rot.length; v++) {

					rotation(w, atom, rot[v]);

					double nb2 = getNumberInter(w);
					double na2 = getNumberAtom(w, dist2);

					if (nb2 < nb) {
						nb = nb2;
						na = na2;
						op1 = v;
					} else {
						if (nb2 == nb) {
							if (na2 > na) {
								nb = nb2;
								na = na2;
								op1 = v;

							}

						}

					}

				}

				rotation(w, atom, (op1 + 1) * 10);

				Bond b = new Bond(atom, w, 1);
				explicit2.add(b);
				water.add(w);

			}

		}

	}

	public void rotation(Atom at, Atom center, double angle) {

		angle = angle * Math.PI / 180;

		double OAx = at.getFx() - center.getFx();
		double OAy = at.getFy() - center.getFy();

		double resx = center.getFx() + Math.cos(angle) * OAx - Math.sin(angle)
				* OAy;
		double resy = center.getFy() + Math.sin(angle) * OAx + Math.cos(angle)
				* OAy;

		at.setFx(resx);
		at.setFy(resy);

	}

	/* creation and traitment of predecesor */
	public Atom traitExplicit(Explicit e, Atom at) {

		Atom res = atom;

		e.setAtLigand(res);

		double dist = 0;

		dist = at.distance(atom);

		Atom ref = at;

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = res.getFx();
		double y1 = res.getFy();
		double x2 = ref.getFx();
		double y2 = ref.getFy();
		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;
		double r = dist;

		if (b != 0.0) {
			double u = a * (c + b * yo) - b * b * xo;
			double v = a * a + b * b;
			double w = c + b * yo;
			double deltap = u * u - v * (b * b * (xo * xo - r * r) + w * w);
			if (deltap >= 0.0) {
				sol1x = (-u + 1.0 * Math.sqrt(deltap)) / v;
				sol1y = -(a * sol1x + c) / b;
				sol2x = (-u + -1.0 * Math.sqrt(deltap)) / v;
				sol2y = -(a * sol2x + c) / b;
			}

		} else {
			sol1x = -c / a;
			double u = sol1x - xo;
			double v = r * r - u * u;
			if (v >= 0.0) {
				sol1y = yo + 1.0 * Math.sqrt(v);
				sol2y = yo + -1.0 * Math.sqrt(v);
				sol2x = sol1x;
			}

		}

		double d1 = Math.sqrt((ref.getFx() - sol1x) * (ref.getFx() - sol1x)
				+ (ref.getFy() - sol1y) * (ref.getFy() - sol1y));
		double d2 = Math.sqrt((ref.getFx() - sol2x) * (ref.getFx() - sol2x)
				+ (ref.getFy() - sol2y) * (ref.getFy() - sol2y));

		double solx = 0;
		double soly = 0;

		if (d1 < d2) {

			solx = sol2x;
			soly = sol2y;

		} else {

			solx = sol1x;
			soly = sol1y;
		}

		double ddx, ddy;

		ddx = solx - at.getFx();
		ddy = soly - at.getFy();

		e.translation(ddx, ddy);

		return res;

	}

	public LinkedList barycenter(LinkedList l) {

		LinkedList res = new LinkedList();

		double sx = 0;
		double sy = 0;

		for (int i = 0; i < l.size(); i++) {
			Atom aa = (Atom) l.get(i);
			sx = sx + aa.getFx();
			sy = sy + aa.getFy();
		}

		sx = sx / l.size();
		sy = sy / l.size();

		res.add(sx);
		res.add(sy);

		return res;
	}

	public LinkedList barycenter2(LinkedList l) {

		LinkedList res = new LinkedList();

		double sx = 0;
		double sy = 0;

		for (int i = 0; i < l.size(); i++) {
			Atom aa = (Atom) l.get(i);
			sx = sx + aa.getX();
			sy = sy + aa.getY();
		}

		sx = sx / l.size();
		sy = sy / l.size();

		res.add(sx);
		res.add(sy);

		return res;
	}

	public boolean resIsInside(Residue r, LinkedList l) {

		Chain c = r.getParent();

		for (int i = 0; i < l.size(); i++) {

			Residue rr = (Residue) l.get(i);
			Chain cc = rr.getParent();

			if ((r.getId() == rr.getId()) && (c.getId() == cc.getId()))
				return true;

		}
		return false;

	}

	public LinkedList getResExpli() {
		return resExplicit;
	}

	public boolean conflit(Explicit e) {

		LinkedList at = e.getAtoms();
		LinkedList bo = e.getBonds();

		for (int j = 0; j < at.size(); j++) {
			Atom a1 = (Atom) at.get(j);
			for (int i = 0; i < atoms.size(); i++) {
				Atom a2 = (Atom) atoms.get(i);

				if (a1.distanceF(a2) < 0.6) {

					return true;
				}
			}
		}

		if (resExplicit.size() > 1) {

			for (int u = 0; u < resExplicit.size(); u++) {

				Explicit ee = (Explicit) resExplicit.get(u);
				if (ee != e) {

					LinkedList tmp = ee.getAtoms();

					for (int j = 0; j < at.size(); j++) {
						Atom a1 = (Atom) at.get(j);
						for (int i = 0; i < tmp.size(); i++) {
							Atom a2 = (Atom) tmp.get(i);
							if (a1.distanceF(a2) < 0.6) {

								return true;
							}
						}
					}
				}
			}

		}

		for (int j = 0; j < bo.size(); j++) {
			Bond b1 = (Bond) bo.get(j);
			for (int i = 0; i < explicit.size(); i++) {
				Bond b2 = (Bond) explicit.get(i);
				if (isIntersect(b1, b2)) {

					return true;
				}
			}
		}

		if (resExplicit.size() != 0) {

			for (int u = 0; u < resExplicit.size(); u++) {

				Explicit ee = (Explicit) resExplicit.get(u);
				if (ee != e) {

					LinkedList tmp = ee.getBonds();

					for (int j = 0; j < bo.size(); j++) {
						Bond b1 = (Bond) bo.get(j);
						for (int i = 0; i < tmp.size(); i++) {
							Bond b2 = (Bond) tmp.get(i);
							if (isIntersect(b1, b2)) {

								return true;
							}
						}
					}
				}
			}

		}

		return false;
	}

	public double getNumberAtom(LigandLight lig, double r) {

		double res = 0;
		int ct = 0;

		LinkedList at = lig.getAtoms();

		for (int j = 0; j < at.size(); j++) {
			Atom a1 = (Atom) at.get(j);
			for (int i = 0; i < atoms.size(); i++) {
				Atom a2 = (Atom) atoms.get(i);

				double dd = a1.distanceF(a2);
				if (dd < r) {

					res = res + dd;
					ct++;
				}
			}
		}

		if (ct == 0)
			return 100;
		else
			res = res / ct;

		return res;
	}

	public double getNumberAtom(Explicit e, double r) {

		double res = 0;
		int ct = 0;

		LinkedList at = e.getAtoms();

		for (int j = 0; j < at.size(); j++) {
			Atom a1 = (Atom) at.get(j);
			for (int i = 0; i < atoms.size(); i++) {
				Atom a2 = (Atom) atoms.get(i);

				double dd = a1.distanceF(a2);
				if (dd < r) {

					res = res + dd;
					ct++;
				}
			}
		}

		for (int u = 0; u < lightLigands.size(); u++) {
			LigandLight lig = (LigandLight) lightLigands.get(u);
			LinkedList att = lig.getAtoms();
			for (int j = 0; j < at.size(); j++) {
				Atom a1 = (Atom) at.get(j);
				for (int i = 0; i < att.size(); i++) {
					Atom a2 = (Atom) att.get(i);

					double dd = a1.distanceF(a2);
					if (dd < r) {

						res = res + dd;
						ct++;
					}
				}
			}
		}

		for (int u = 0; u < resExplicit.size(); u++) {
			Explicit ee = (Explicit) resExplicit.get(u);
			LinkedList att = ee.getAtoms();
			for (int j = 0; j < at.size(); j++) {
				Atom a1 = (Atom) at.get(j);
				for (int i = 0; i < att.size(); i++) {
					Atom a2 = (Atom) att.get(i);

					double dd = a1.distanceF(a2);
					if (dd < r) {

						res = res + dd;
						ct++;
					}
				}
			}
		}

		if (ct == 0)
			return 100;
		else
			res = res / ct;

		return res;
	}

	public double getNumberAtom(Atom w, double r) {

		double res = 0;
		int ct = 0;

		for (int i = 0; i < atoms.size(); i++) {
			Atom a2 = (Atom) atoms.get(i);

			double dd = w.distanceF(a2);
			if (dd < r) {

				res = res + dd;
				ct++;
			}
		}

		for (int u = 0; u < lightLigands.size(); u++) {
			LigandLight lig = (LigandLight) lightLigands.get(u);
			LinkedList att = lig.getAtoms();

			for (int i = 0; i < att.size(); i++) {
				Atom a2 = (Atom) att.get(i);

				double dd = w.distanceF(a2);
				if (dd < r) {

					res = res + dd;
					ct++;
				}

			}
		}

		for (int u = 0; u < resExplicit.size(); u++) {
			Explicit ee = (Explicit) resExplicit.get(u);
			LinkedList att = ee.getAtoms();

			for (int i = 0; i < att.size(); i++) {
				Atom a2 = (Atom) att.get(i);

				double dd = w.distanceF(a2);
				if (dd < r) {

					res = res + dd;
					ct++;
				}

			}
		}

		for (int u = 0; u < water.size(); u++) {
			Atom w2 = (Atom) water.get(u);

			if (w.getId() != w2.getId()) {

				double dd = w.distanceF(w2);
				if (dd < r) {

					res = res + dd;
					ct++;
				}
			}
		}

		if (ct == 0)
			return 100;
		else
			res = res / ct;

		return res;
	}

	public double getNumberInter(LigandLight lig) {

		LinkedList bo = lig.getBonds();

		double res = 0;

		for (int j = 0; j < bo.size(); j++) {
			Bond b1 = (Bond) bo.get(j);
			for (int i = 0; i < explicit.size(); i++) {
				Bond b2 = (Bond) explicit.get(i);
				Atom a1 = b2.getFirst();
				Atom a2 = b2.getSecond();
				if ((a1.getFx() != 0) && (a2.getFx() != 0)) {
					if (isIntersect(b1, b2)) {
						res = res + 1;
					}
				}
			}
		}

		return res;

	}

	public double getNumberInter(Explicit e) {

		LinkedList bo = e.getBonds();

		double res = 0;

		for (int j = 0; j < bo.size(); j++) {
			Bond b1 = (Bond) bo.get(j);
			for (int i = 0; i < explicit.size(); i++) {
				Bond b2 = (Bond) explicit.get(i);
				Atom a1 = b2.getFirst();
				Atom a2 = b2.getSecond();
				if ((a1.getFx() != 0) && (a2.getFx() != 0)) {
					if (isIntersect(b1, b2)) {
						res = res + 1;
					}
				}
			}
		}

		for (int u = 0; u < lightLigands.size(); u++) {
			LigandLight lig = (LigandLight) lightLigands.get(u);
			LinkedList bb = lig.getBonds();

			for (int j = 0; j < bo.size(); j++) {
				Bond b1 = (Bond) bo.get(j);
				for (int i = 0; i < bb.size(); i++) {
					Bond b2 = (Bond) bb.get(i);

					if (isIntersect(b1, b2)) {
						res = res + 1;
					}

				}
			}

		}

		for (int u = 0; u < resExplicit.size(); u++) {
			Explicit ee = (Explicit) resExplicit.get(u);
			LinkedList bb = ee.getBonds();

			for (int j = 0; j < bo.size(); j++) {
				Bond b1 = (Bond) bo.get(j);
				for (int i = 0; i < bb.size(); i++) {
					Bond b2 = (Bond) bb.get(i);

					if (isIntersect(b1, b2)) {
						res = res + 1;
					}

				}
			}

		}

		return res;

	}

	public double getNumberInter(Atom w) {

		Bond b1 = new Bond(atom, w, 1);

		double res = 0;

		for (int i = 0; i < explicit.size(); i++) {
			Bond b2 = (Bond) explicit.get(i);
			Atom a1 = b2.getFirst();
			Atom a2 = b2.getSecond();
			if ((a1.getFx() != 0) && (a2.getFx() != 0)) {
				if (isIntersect(b1, b2)) {
					res = res + 1;
				}
			}

		}

		for (int u = 0; u < lightLigands.size(); u++) {
			LigandLight lig = (LigandLight) lightLigands.get(u);
			LinkedList bb = lig.getBonds();

			for (int i = 0; i < bb.size(); i++) {
				Bond b2 = (Bond) bb.get(i);

				if (isIntersect(b1, b2)) {
					res = res + 1;
				}

			}

		}

		for (int u = 0; u < resExplicit.size(); u++) {
			Explicit ee = (Explicit) resExplicit.get(u);
			LinkedList bb = ee.getBonds();

			for (int i = 0; i < bb.size(); i++) {
				Bond b2 = (Bond) bb.get(i);

				if (isIntersect(b1, b2)) {
					res = res + 1;
				}

			}

		}

		return res;

	}

	public boolean isIntersect(Bond bond1, Bond bond2) {

		Atom aa1 = bond1.getFirst();
		Atom aa2 = bond1.getSecond();

		Atom bb1 = bond2.getFirst();
		Atom bb2 = bond2.getSecond();

		double Ax = aa1.getFx();
		double Ay = aa1.getFy();
		double Bx = aa2.getFx();
		double By = aa2.getFy();
		double Cx = bb1.getFx();
		double Cy = bb1.getFy();
		double Dx = bb2.getFx();
		double Dy = bb2.getFy();

		double Sx;
		double Sy;

		if (Ax == Bx) {
			if (Cx == Dx)
				return false;
			else {
				double pCD = (Cy - Dy) / (Cx - Dx);
				Sx = Ax;
				Sy = pCD * (Ax - Cx) + Cy;
			}
		} else {
			if (Cx == Dx) {
				double pAB = (Ay - By) / (Ax - Bx);
				Sx = Cx;
				Sy = pAB * (Cx - Ax) + Ay;
			} else {
				double pCD = (Cy - Dy) / (Cx - Dx);
				double pAB = (Ay - By) / (Ax - Bx);
				double oCD = Cy - pCD * Cx;
				double oAB = Ay - pAB * Ax;
				Sx = (oAB - oCD) / (pCD - pAB);
				Sy = pCD * Sx + oCD;
			}
		}
		if ((Sx <= Ax && Sx <= Bx) | (Sx >= Ax && Sx >= Bx)
				| (Sx <= Cx && Sx <= Dx) | (Sx >= Cx && Sx >= Dx)
				| (Sy <= Ay && Sy <= By) | (Sy >= Ay && Sy >= By)
				| (Sy <= Cy && Sy <= Dy) | (Sy >= Cy && Sy >= Dy))
			return false;
		return true;

	}

}
