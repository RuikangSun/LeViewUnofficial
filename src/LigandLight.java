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
 This class defines just a ligand (without contacts, HBonds...).
 */

import java.io.*;
import java.util.*;

public class LigandLight {

	Chain chain;
	LinkedList residues;
	LinkedList atoms;
	LinkedList bonds;
	LinkedList explicit;
	/* Hbonds contains the list of Hbond */
	LinkedList Hbonds;
	/* contacts contains the list of residues closed to the ligand */
	LinkedList contacts;
	/* close contains all the residues that have an atom close to the ligand */
	LinkedList close;
	LinkedList copy;
	int cutoff;
	Molecule mol;
	private LinkedList listTemp;
	private LinkedList forbidden;
	LinkedList chaines;
	LinkedList forbiddenEx;
	LinkedList contactList;
	LinkedList water;
	LinkedList HbondsW;
	LinkedList posName;
	LinkedList rings;
	LinkedList resExplicit;

	boolean selected;

	public LigandLight() {
	}

	public LigandLight(Chain chain) {

		this.chain = chain;
		mol = chain.getParent();
		cutoff = 4;
		residues = chain.getResidues();
		atoms = new LinkedList();
		bonds = new LinkedList();
		rings = new LinkedList();
		forbidden = new LinkedList();
		chaines = new LinkedList();
		water = mol.getWater();
		HbondsW = new LinkedList();
		explicit = new LinkedList();
		forbiddenEx = new LinkedList();
		contactList = new LinkedList();
		posName = new LinkedList();

		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);
			LinkedList at = res.getAtoms();
			for (int j = 0; j < at.size(); j++) {
				Atom atom = (Atom) at.get(j);
				atom.init();

				LinkedList l1 = atom.getNeig();

				for (int v = 0; v < l1.size(); v++) {
					Atom a = (Atom) l1.get(v);
					a.init();
				}

				LinkedList l2 = atom.getHneig();

				for (int v = 0; v < l2.size(); v++) {
					Atom a = (Atom) l2.get(v);
					a.init();
				}

				atoms.add(atom);
			}
		}

		updateConect();

		/* create atom and bond lists */
		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);
			LinkedList at = res.getAtoms();
			for (int j = 0; j < at.size(); j++) {
				Atom atom = (Atom) at.get(j);

				LinkedList neig = atom.getNeig();
				LinkedList Hneig = atom.getHneig();

				for (int u = 0; u < neig.size(); u++) {
					Atom att = (Atom) neig.get(u);

					Bond b = new Bond(atom, att, 1);
					if (!bondExist(b)) {

						bonds.add(b);
					}

				}

			}
		}
		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);
			LinkedList at = res.getAtoms();
			for (int j = 0; j < at.size(); j++) {
				Atom atom = (Atom) at.get(j);
				LinkedList neig = atom.getNeig();
				LinkedList Hneig = atom.getHneig();

			}

		}

		updateBonds();

		setTerm();

		getRings();
		updateConjugate();
		updateAromatic();

		getChaines();
		updateTypeAtom();
		updateCoord();

		for (int i = 0; i < chaines.size(); i++) {
			Chaine a = (Chaine) chaines.get(i);

		}

		layout();

		updateAromBond();

		resolutionConflit();

		posName = getPoseName();

	}

	public LinkedList getAtoms() {
		return atoms;
	}

	public void translation(double x, double y) {

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			double ddx, ddy;

			ddx = at.getFx() + x;
			ddy = at.getFy() + y;

			at.setFx(ddx);
			at.setFy(ddy);

		}

	}

	public void setSelected(boolean b) {
		selected = b;
	}

	public boolean getSelected() {
		return selected;
	}

	public void rotation(Atom centre, double angle) {

		angle = Math.PI * angle / 180;

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);

			double OAx = at.getFx() - centre.getFx();
			double OAy = at.getFy() - centre.getFy();

			double resx = centre.getFx() + Math.cos(angle) * OAx
					- Math.sin(angle) * OAy;
			double resy = centre.getFy() + Math.sin(angle) * OAx
					+ Math.cos(angle) * OAy;

			at.setFx(resx);
			at.setFy(resy);

		}

	}

	public LigandLight getCopy() {
		LigandLight ex = new LigandLight();
		ex.setResidues(residues);
		ex.setParent(chain);
		ex.setPosName(posName);
		LinkedList temp = new LinkedList();
		/* copying atoms */
		for (int i = 0; i < atoms.size(); i++) {
			Atom a1 = (Atom) atoms.get(i);
			temp.add(a1.getCopy());
		}
		ex.setAtoms(temp);
		temp = new LinkedList();
		/* copying bonds */
		for (int i = 0; i < bonds.size(); i++) {
			Bond b = (Bond) bonds.get(i);
			if ((getCorres(b.getFirst(), ex.getAtoms()) != null)
					&& (getCorres(b.getSecond(), ex.getAtoms()) != null)) {
				Bond b1 = new Bond(getCorres(b.getFirst(), ex.getAtoms()),
						getCorres(b.getSecond(), ex.getAtoms()), b.getOrder());
				temp.add(b1);
			}
		}
		ex.setBonds(temp);
		ex.setRing(rings);

		return ex;

	}

	public void setRing(LinkedList l) {
		rings = l;
	}

	public Atom getCorres(Atom a, LinkedList l) {
		for (int i = 0; i < l.size(); i++) {
			Atom at = (Atom) l.get(i);
			if (at.getId() == a.getId())
				return at;
		}
		return null;
	}

	public void setAtoms(LinkedList l) {
		atoms = l;
	}

	public void setBonds(LinkedList l) {
		bonds = l;
	}

	public void setPosName(LinkedList l) {
		posName = l;
	}

	public LinkedList getRing() {
		return rings;
	}

	public Molecule getMol() {
		return mol;
	}

	public Chain getParent() {
		return chain;
	}

	public void setParent(Chain c) {
		chain = c;
	}

	public void setResidues(LinkedList l) {
		residues = l;
	}

	public LinkedList getContacts() {
		return contactList;
	}

	void updateConect() {

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			LinkedList neig = at.getNeig();

			for (int j = 0; j < neig.size(); j++) {
				Atom a = (Atom) neig.get(j);
				if (getAtom(a.getId()) == null) {

					at.removeNeig(a.getId());
					at.addHneig(a);
				}

			}

		}

	}

	public LinkedList getPoseName() {

		LinkedList res = new LinkedList();

		for (int i = 0; i < residues.size(); i++) {
			Residue r = (Residue) residues.get(i);

			res.add(getPosName(r));

		}

		return res;

	}

	public LinkedList returnPoseName() {
		return posName;
	}

	public LinkedList getHbond() {
		return Hbonds;
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

		/* intersection bonds */
		for (int j = 0; j < bo.size(); j++) {
			Bond b1 = (Bond) bo.get(j);
			for (int i = 0; i < bonds.size(); i++) {
				Bond b2 = (Bond) bonds.get(i);
				if (isIntersect(b1, b2)) {

					return true;
				}
			}
		}

		if (resExplicit.size() > 1) {

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

	public void resolutionConflit() {

		if (!getConflit())
			return;

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			if (at.getPivot()) {

				LinkedList atomPivot = getListRota(at);
				Rotation(at, atomPivot, 120);
				if (getConflit()) {
					Rotation(at, atomPivot, 120);
					if (getConflit()) {
						Rotation(at, atomPivot, 120);
					} else {
						return;
					}
				} else {
					return;
				}

			}

		}

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);

			if ((at.getType() == 2) && (at.getNeig().size() == 3)) {

				LinkedList atomPivot = getListRota2(at);
				LinkedList l1 = (LinkedList) atomPivot.get(0);
				LinkedList l2 = (LinkedList) atomPivot.get(1);
				Rotation(at, l1, 120);
				Rotation(at, l2, -120);
				if (getConflit()) {
					Rotation(at, l1, -240);
					Rotation(at, l2, 240);
					if (getConflit()) {
						Rotation(at, l1, 120);
						Rotation(at, l2, -120);
					} else {
						return;
					}
				} else {
					return;
				}

			}

		}

		for (int i = 0; i < rings.size(); i++) {
			Ring r = (Ring) rings.get(i);
			if (r.getConjugate().size() == 0) {
				LinkedList l = r.getAtoms();
				Atom d = null;
				for (int j = 0; j < l.size(); j++) {
					Atom t = getAtom((Integer) l.get(j));
					if (t.getNeig().size() > 2) {
						d = t;
						break;
					}
				}
				if (d != null) {
					symetrieRing(d, r);
					if (getConflit()) {
						symetrieRing(d, r);
					} else {
						return;
					}
				}
			}
		}

	}

	public void symetrieRing(Atom dep, Ring ring) {

		int n = ring.getLenght();
		Atom[] ringAtoms = new Atom[n];
		LinkedList temp = ring.getAtoms();

		int pos = 0;

		for (int i = 0; i < temp.size(); i++) {
			if ((Integer) temp.get(i) == dep.getId()) {
				pos = i;
				break;
			}
		}

		ringAtoms[0] = getAtom((Integer) temp.get(pos));

		int cont = 1;

		for (int i = pos + 1; i < temp.size(); i++) {
			ringAtoms[cont] = getAtom((Integer) temp.get(i));
			cont++;
		}
		for (int i = 0; i < pos; i++) {
			ringAtoms[cont] = getAtom((Integer) temp.get(i));
			cont++;
		}

		Atom d1;
		Atom d2;

		if (n % 2 == 0) {
			d1 = ringAtoms[0];
			d2 = ringAtoms[n / 2];
		} else {
			d1 = ringAtoms[0];
			d2 = d1.clone();
			double tx = (ringAtoms[n / 2].getFx() + ringAtoms[n / 2 + 1]
					.getFx()) / 2;
			double ty = (ringAtoms[n / 2].getFy() + ringAtoms[n / 2 + 1]
					.getFy()) / 2;
			d2.setFx(tx);
			d2.setFy(ty);
		}

		double a = d2.getFy() - d1.getFy();
		double b = d1.getFx() - d2.getFx();
		double c = d1.getFy() * d2.getFx() - d1.getFx() * d2.getFy();

		for (int i = 0; i < n; i++) {
			Atom at = ringAtoms[i];

			if ((at.getNeig().size() > 2) && (i != 0)) {
				LinkedList tmp = getListRota2(ringAtoms[i]);
				LinkedList l1 = (LinkedList) tmp.get(0);
				LinkedList l2 = (LinkedList) tmp.get(1);

				LinkedList lll = new LinkedList();

				if (i != n - 1) {
					if ((!isInside(ringAtoms[i + 1], l1))
							&& (!isInside(ringAtoms[i - 1], l1)))
						lll = l1;
					else
						lll = l2;
				} else {
					if ((!isInside(ringAtoms[0], l1))
							&& (!isInside(ringAtoms[i - 1], l1)))
						lll = l1;
					else
						lll = l2;
				}

				for (int u = 0; u < lll.size(); u++) {
					Atom tt = (Atom) lll.get(u);

					sym(tt, a, b, c);
				}
			}

			sym(at, a, b, c);

		}

	}

	public void sym(Atom at, double a, double b, double c) {

		double Ax = at.getFx();
		double Ay = at.getFy();

		double u = (a * at.getFx() + b * at.getFy() + c) / (a * a + b * b);
		at.setFx(at.getFx() - u * a);
		at.setFy(at.getFy() - u * b);

		at.setFx(2.0 * at.getFx() - Ax);
		at.setFy(2.0 * at.getFy() - Ay);

	}

	public boolean getConflit() {

		for (int i = 0; i < atoms.size() - 1; i++) {
			Atom a1 = (Atom) atoms.get(i);

			for (int j = i + 1; j < atoms.size(); j++) {
				Atom a2 = (Atom) atoms.get(j);

				if (a1.distanceF(a2) < 0.6)
					return true;

			}

		}
		return false;
	}

	public boolean conflit(Hbond h, LinkedList l) {

		Atom atom = h.getAtomRes();

		for (int i = 0; i < atoms.size(); i++) {
			Atom a2 = (Atom) atoms.get(i);

			if (atom.distanceF(a2) < 1.2) {

				return true;
			}
		}

		if (resExplicit.size() > 0) {

			for (int u = 0; u < resExplicit.size(); u++) {

				Explicit ee = (Explicit) resExplicit.get(u);

				LinkedList tmp = ee.getAtoms();

				for (int i = 0; i < tmp.size(); i++) {
					Atom a2 = (Atom) tmp.get(i);
					if (atom.distanceF(a2) < 0.6) {

						return true;
					}
				}
			}

		}

		for (int i = 0; i < l.size(); i++) {
			Hbond hb = (Hbond) l.get(i);
			Atom atom2 = hb.getAtomRes();

			if (h.equals(hb))
				break;

			if (atom.distanceF(atom2) < 0.6) {

				return true;
			}

		}

		/* intersection bonds */

		for (int i = 0; i < bonds.size(); i++) {
			Bond b2 = (Bond) bonds.get(i);
			if (isIntersect(h, b2)) {

				return true;
			}
		}

		if (resExplicit.size() > 0) {

			for (int u = 0; u < resExplicit.size(); u++) {

				Explicit ee = (Explicit) resExplicit.get(u);

				LinkedList tmp = ee.getBonds();

				for (int i = 0; i < tmp.size(); i++) {
					Bond b2 = (Bond) tmp.get(i);
					if (isIntersect(h, b2)) {

						return true;
					}
				}
			}

		}

		for (int i = 0; i < l.size(); i++) {
			Hbond hb = (Hbond) l.get(i);
			Atom atom2 = hb.getAtomRes();

			if (h.equals(hb))
				break;

			if (isIntersect(h, hb)) {

				return true;
			}

		}

		return false;
	}

	public LinkedList getContactList() {
		LinkedList res = new LinkedList();

		for (int i = 0; i < contacts.size(); i++) {
			Residue r = (Residue) contacts.get(i);

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
			LinkedList lll = ref.getNeig();
			Atom res = (Atom) lll.get(0);

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

		for (int i = 0; i < Hbonds.size(); i++) {
			Hbond hb = (Hbond) Hbonds.get(i);
			Atom a2 = hb.getAtomRes();

			double dd = Math.sqrt((atom.getFx() - a2.getFx())
					* (atom.getFx() - a2.getFx()) + (atom.getFy() - a2.getFy())
					* (atom.getFy() - a2.getFy()));

			if (dd <= r) {

				return true;
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

	public int getNumberInter(Hbond h, LinkedList l) {

		int res = 0;

		Atom atom = h.getAtomRes();

		/* intersection bonds */

		for (int i = 0; i < bonds.size(); i++) {
			Bond b2 = (Bond) bonds.get(i);
			if (isIntersect(h, b2)) {

				res++;
			}
		}

		if (resExplicit.size() > 0) {

			for (int u = 0; u < resExplicit.size(); u++) {

				Explicit ee = (Explicit) resExplicit.get(u);

				LinkedList tmp = ee.getBonds();

				for (int i = 0; i < tmp.size(); i++) {
					Bond b2 = (Bond) tmp.get(i);
					if (isIntersect(h, b2)) {

						res++;
					}
				}
			}

		}

		for (int i = 0; i < explicit.size(); i++) {

			Bond b2 = (Bond) explicit.get(i);

			if (isIntersect(h, b2)) {

				res++;
			}

		}

		for (int i = 0; i < l.size(); i++) {
			Hbond hb = (Hbond) l.get(i);
			Atom atom2 = hb.getAtomRes();

			if (h.equals(hb))
				break;

			if (isIntersect(h, hb)) {

				res = res + 1;
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

	public boolean isIntersect(Hbond bond1, Bond bond2) {

		Atom aa1 = bond1.getDon();
		Atom aa2 = bond1.getAcc();

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

	public boolean isIntersect(Hbond bond1, Hbond bond2) {

		Atom aa1 = bond1.getDon();
		Atom aa2 = bond1.getAcc();

		Atom bb1 = bond2.getDon();
		Atom bb2 = bond2.getAcc();

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

	public void getResExplicit() {

		LinkedList res = new LinkedList();
		LinkedList aaa = new LinkedList();

		for (int i = 0; i < explicit.size(); i++) {
			Bond b = (Bond) explicit.get(i);
			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();

			Residue r1 = a1.getParent();
			Residue r2 = a2.getParent();

			if ((!resIsInside(r1, residues)) && (!resIsInside(r1, res))) {
				res.add(r1);
				aaa.add(a1);
			}
			if ((!resIsInside(r2, residues)) && (!resIsInside(r2, res))) {
				res.add(r2);
				aaa.add(a2);
			}

		}

		for (int i = 0; i < res.size(); i++) {
			Atom att = (Atom) aaa.get(i);
			Residue rrr = (Residue) res.get(i);

			Explicit e = new Explicit(rrr, att);
			forbiddenEx.add(rrr);
			traitExplicit(e, att);

			double[] rot = new double[36];
			for (int s = 0; s < rot.length; s++) {
				rot[s] = 10;
			}

			boolean b = false;

			for (int v = 0; v < rot.length; v++) {

				if (conflit(e)) {

					e.rotation(e.getAtLigand(), rot[v]);
				} else {
					b = true;

					break;
				}

			}

			if (!b) {

				for (int v = 0; v < rot.length; v++) {

					if (conflit(e)) {

						e.rotation(att, rot[v]);

						for (int u = 0; u < rot.length; u++) {

							if (conflit(e)) {

								e.rotation(e.getAtLigand(), rot[u]);
							} else {
								b = true;

								break;
							}

						}

					} else {
						b = true;

						break;
					}

				}

			}

			resExplicit.add(e);
		}

	}

	/* creation of predecesor atom */
	public Atom traitExplicit(Explicit e, Atom at) {

		Atom res = null;
		res = at.clone();
		LinkedList voisin = new LinkedList();

		for (int i = 0; i < atoms.size(); i++) {
			Atom aa = (Atom) atoms.get(i);
			LinkedList HN = aa.getHneig();

			for (int j = 0; j < HN.size(); j++) {
				Atom aaa = (Atom) HN.get(j);
				if (aaa.getId() == at.getId())
					voisin.add(aa);
			}
		}

		if (voisin.size() == 0)
			return res;

		if (voisin.size() == 1) {
			res = (Atom) voisin.get(0);
		} else {

			LinkedList tmp = barycenter(voisin);
			res.setFx((Double) tmp.get(0));
			res.setFy((Double) tmp.get(1));
		}

		e.setAtLigand(res);

		double dist = 0;

		for (int i = 0; i < voisin.size(); i++) {
			Atom aaa = (Atom) voisin.get(i);

			dist = dist + at.distance(aaa);

		}

		dist = dist / voisin.size();

		Atom temp = (Atom) voisin.get(0);

		LinkedList l = temp.getNeig();
		Atom ref = null;
		if (l.size() == 0) {
			ref = at;
		} else {

			ref = (Atom) l.get(0);
		}

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

	public void updateTypeAtom() {
		for (int i = 0; i < rings.size(); i++) {
			Ring r = (Ring) rings.get(i);
			int type;
			if (r.isAromatic())
				type = 4;
			else
				type = 3;
			LinkedList att = r.getAtoms();
			for (int j = 0; j < att.size(); j++) {
				int id = (Integer) att.get(j);
				Atom at = getAtom(id);
				at.setType(type);
			}
		}

		for (int i = 0; i < chaines.size(); i++) {
			Chaine r = (Chaine) chaines.get(i);

			LinkedList att = r.getAtoms();
			for (int j = 0; j < att.size(); j++) {
				int id = (Integer) att.get(j);
				Atom at = getAtom(id);
				at.setType(1);
			}
		}

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			if (at.getType() == -1) {
				if (at.getTerm())
					at.setType(6);
				else
					at.setType(5);
			}
		}

	}

	public void updateCoord() {
		for (int i = 0; i < atoms.size(); i++) {
			Atom atom = (Atom) atoms.get(i);
			atom.setFx(atom.getX());
			atom.setFy(atom.getY());
		}
	}

	public void getRings() {
		for (int i = 0; i < atoms.size(); i++) {
			Atom atom = (Atom) atoms.get(i);
			if (!estDedans(atom.getId(), forbidden)) {
				LinkedList l = lance(atom.getId());
				if (l.size() != 0) {

					for (int u = 0; u < l.size(); u++) {
						int id = (Integer) l.get(u);

						forbidden.add(id);
					}
					Ring r = new Ring(l);
					rings.add(r);
				}
			}
		}

	}

	public void updateConjugate() {
		for (int i = 0; i < rings.size() - 1; i++) {
			Ring r1 = (Ring) rings.get(i);
			for (int j = i + 1; j < rings.size(); j++) {
				Ring r2 = (Ring) rings.get(j);
				if (r1.isInside(r2.getAtoms())) {
					r1.addConjugate(r2);
					r2.addConjugate(r1);
				}
			}
		}
	}

	public LinkedList getAtomLinkedChain(Chaine c) {
		LinkedList l = new LinkedList();
		LinkedList cc = c.getAtoms();

		for (int i = 0; i < cc.size(); i++) {
			int id = (Integer) cc.get(i);
			Atom at = getAtom(id);
			LinkedList neig = at.getNeig();

			for (int j = 0; j < neig.size(); j++) {
				Atom a = (Atom) neig.get(j);
				int ID = a.getId();
				if (cc.indexOf(ID) == -1)
					l.add(a);
			}
		}

		return l;

	}

	public void setTerm() {

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);

			if (at.getNeig().size() == 1)
				at.setTerm(true);
		}
	}

	public LinkedList getResidues() {
		return residues;
	}

	public LinkedList getBonds() {
		return bonds;
	}

	public LinkedList getExplicit() {
		return explicit;
	}

	public LinkedList getClose() {
		LinkedList l = new LinkedList();
		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);
			LinkedList ll = getCloseRes(res);
			for (int h = 0; h < ll.size(); h++) {
				Residue rrr = (Residue) ll.get(h);
				if ((!isInside(rrr, l)) && (!isInside(rrr, forbiddenEx))) {
					l.add(rrr);
				}
			}
		}
		return l;
	}

	public LinkedList lance(int u) {
		LinkedList res = new LinkedList();
		listTemp = new LinkedList();
		LinkedList open = new LinkedList();
		LinkedList closed = new LinkedList();
		closed.add(u);
		Atom at = getAtom(u);

		for (int i = 0; i < at.getNeig().size(); i++) {
			Atom atom = (Atom) at.getNeig().get(i);
			open.add(atom.getId());
		}

		trait(closed, open, u);

		int c = 100;

		for (int i = 0; i < listTemp.size(); i++) {
			LinkedList l = (LinkedList) listTemp.get(i);
			if (l.size() < c) {
				res = l;
				c = l.size();
			}
		}

		return res;
	}

	public void trait(LinkedList closed, LinkedList open, int pos) {

		LinkedList temp = new LinkedList();

		if (open.size() == 0) {
		} else {
			do {

				int element = (Integer) open.get(0);
				open.remove(0);
				temp.clear();
				Atom at = getAtom(element);

				for (int i = 0; i < at.getNeig().size(); i++) {

					Atom atom = (Atom) at.getNeig().get(i);

					temp.add(atom.getId());
				}

				if ((estDedans(pos, temp)) && (closed.size() > 3)) {
					LinkedList res = new LinkedList();
					closed.add(element);
					for (int s = 0; s < closed.size(); s++) {
						res.add((Integer) closed.get(s));
					}

					listTemp.add(res);
				} else {
					for (int i = 0; i < temp.size(); i++) {

						if (estDedans((Integer) temp.get(i), closed)) {
							temp.remove(i);
							i--;
						}
					}

					closed.add(element);

					trait(closed, temp, pos);

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

	LinkedList resChain = new LinkedList();
	int max = 0;

	public void getChaines() {
		LinkedList open = new LinkedList();
		for (int i = 0; i < atoms.size(); i++) {
			Atom atom = (Atom) atoms.get(i);
			if ((atom.getElement().equals("C"))
					|| (atom.getElement().equals("N"))
					|| (atom.getElement().equals("O"))
					|| (atom.getElement().equals("S"))) {
				if (!atom.getTerm()) {
					if (!estDedans(atom.getId(), forbidden)) {
						open.add(atom.getId());
					}
				}
			}
		}
		if (open.size() != 0) {
			do {
				lance2(open);
				open = new LinkedList();
				for (int i = 0; i < atoms.size(); i++) {
					Atom atom = (Atom) atoms.get(i);
					if ((atom.getElement().equals("C"))
							|| (atom.getElement().equals("N"))
							|| (atom.getElement().equals("O"))
							|| (atom.getElement().equals("S"))) {
						if (!atom.getTerm()) {
							if (!estDedans(atom.getId(), forbidden)) {
								open.add(atom.getId());
							}
						}
					}
				}
			} while (open.size() != 0);
		}
	}

	LinkedList resChainBis;

	public void lance2(LinkedList open) {
		resChainBis = new LinkedList();
		resChain = new LinkedList();
		max = 0;

		LinkedList open2 = new LinkedList();
		for (int i = 0; i < open.size(); i++) {
			open2.add((Integer) open.get(i));
		}

		LinkedList closed = new LinkedList();

		closed = new LinkedList();
		trait2(closed, open);

		if (resChain.size() != 0) {
			if (resChain.size() == 1) {
				int e = (Integer) resChain.get(0);
				Atom at = getAtom(e);
				at.setType(5);
				forbidden.add(e);
			} else {
				if (resChain.size() < 3) {
					for (int i = 0; i < resChain.size(); i++) {
						int e = (Integer) resChain.get(i);
						Atom at = getAtom(e);
						at.setType(2);
						forbidden.add(e);
					}
				} else {
					trait2Bis(closed, open2, resChain.size());

					if (resChainBis.size() == 0) {
						for (int i = 0; i < resChain.size(); i++) {
							int e = (Integer) resChain.get(i);

							forbidden.add(e);
						}
						Chaine chaine = new Chaine(resChain);
						chaines.add(chaine);
					} else {

						int resPos = 0;
						int cnt = 0;

						for (int i = 0; i < resChainBis.size(); i++) {
							LinkedList rct = (LinkedList) resChainBis.get(i);

							int cnt2 = 0;
							for (int j = 0; j < rct.size(); j++) {
								Atom att = getAtom((Integer) rct.get(j));
								String natt = att.getName();

								if ((natt.equals("C")) || (natt.equals("CA"))
										|| (natt.equals("O"))
										|| (natt.equals("N"))) {
									cnt2 = cnt2 + 1;
								}
								if (cnt2 > cnt) {
									cnt = cnt2;
									resPos = i;
								}

							}

						}

						LinkedList resUlt = (LinkedList) resChainBis
								.get(resPos);

						for (int i = 0; i < resUlt.size(); i++) {
							int e = (Integer) resUlt.get(i);

							forbidden.add(e);
						}
						Chaine chaine = new Chaine(resUlt);
						chaines.add(chaine);

					}

				}
			}
		}

	}

	public void trait2(LinkedList closed, LinkedList open) {

		LinkedList temp = new LinkedList();

		if (open.size() == 0) {
		} else {
			do {

				int element = (Integer) open.get(0);
				open.remove(0);
				temp.clear();
				closed.add(element);
				Atom at = getAtom(element);
				for (int i = 0; i < at.getNeig().size(); i++) {
					Atom atom = (Atom) at.getNeig().get(i);
					if ((atom.getElement().equals("C"))
							|| (atom.getElement().equals("N"))
							|| (atom.getElement().equals("O"))
							|| (atom.getElement().equals("S"))) {
						if (!atom.getTerm()) {
							if (!estDedans(atom.getId(), forbidden)) {
								temp.add(atom.getId());
							}
						}
					}
				}

				if (closed.size() > max) {

					max = closed.size();
					resChain = new LinkedList();

					for (int s = 0; s < closed.size(); s++) {
						resChain.add((Integer) closed.get(s));
					}

				}
				for (int i = 0; i < temp.size(); i++) {

					if (estDedans((Integer) temp.get(i), closed)) {
						temp.remove(i);
						i--;
					}
				}

				trait2(closed, temp);

				closed.remove(closed.size() - 1);
			} while (open.size() != 0);

		}

	}

	public void trait2Bis(LinkedList closed, LinkedList open, int size) {

		LinkedList temp = new LinkedList();

		if (open.size() == 0) {
		} else {
			do {

				int element = (Integer) open.get(0);
				open.remove(0);
				temp.clear();
				closed.add(element);
				Atom at = getAtom(element);
				for (int i = 0; i < at.getNeig().size(); i++) {
					Atom atom = (Atom) at.getNeig().get(i);
					if ((atom.getElement().equals("C"))
							|| (atom.getElement().equals("N"))
							|| (atom.getElement().equals("O"))
							|| (atom.getElement().equals("S"))) {
						if (!atom.getTerm()) {
							if (!estDedans(atom.getId(), forbidden)) {
								temp.add(atom.getId());
							}
						}
					}
				}

				if (closed.size() == size) {

					resChain = new LinkedList();

					for (int s = 0; s < closed.size(); s++) {
						resChain.add((Integer) closed.get(s));
					}

					resChainBis.add(resChain);

				}
				for (int i = 0; i < temp.size(); i++) {

					if (estDedans((Integer) temp.get(i), closed)) {
						temp.remove(i);
						i--;
					}
				}

				trait2Bis(closed, temp, size);

				closed.remove(closed.size() - 1);
			} while (open.size() != 0);

		}

	}

	public void layout() {
		Atom first = getFirst();

		LinkedList neig = first.getNeig();
		LinkedList close = new LinkedList();
		LinkedList open = new LinkedList();
		for (int i = 0; i < neig.size(); i++) {
			Atom at = (Atom) neig.get(i);
			open.add(at);
		}
		close.add(first);
		if (!first.getFlag()) {

			if (first.getType() == 1) {
				traitChaine(first, getChaine(first));
			}
			if ((first.getType() == 2) || (first.getType() == 5)
					|| (first.getType() == 6))
				traitAtom(first);
			if ((first.getType() == 3) || (first.getType() == 4))
				traitRing(first, getRing(first));
		}

		traitLayout(close, open);
	}

	public void traitLayout(LinkedList closed, LinkedList open) {

		LinkedList temp = new LinkedList();

		if (open.size() == 0) {
		} else {
			do {

				Atom at = (Atom) open.get(0);
				open.remove(0);
				temp.clear();
				closed.add(at);

				if (!at.getFlag()) {

					if (at.getType() == 1) {
						traitChaine(at, getChaine(at));
					}
					if ((at.getType() == 2) || (at.getType() == 5)
							|| (at.getType() == 6))
						traitAtom(at);
					if ((at.getType() == 3) || (at.getType() == 4))
						traitRing(at, getRing(at));
				}

				for (int i = 0; i < at.getNeig().size(); i++) {
					Atom atom = (Atom) at.getNeig().get(i);
					temp.add(atom);
				}

				for (int i = 0; i < temp.size(); i++) {

					if (estDedans2((Atom) temp.get(i), closed)) {
						temp.remove(i);
						i--;
					}
				}

				traitLayout(closed, temp);

			} while (open.size() != 0);

		}

	}

	public boolean estDedans2(Atom at, LinkedList l) {
		for (int i = 0; i < l.size(); i++) {
			Atom atom = (Atom) l.get(i);
			if (at.getId() == atom.getId())
				return true;
		}
		return false;
	}

	public Chaine getChaine(Atom at) {
		for (int i = 0; i < chaines.size(); i++) {
			Chaine c = (Chaine) chaines.get(i);
			LinkedList atoms = c.getAtoms();
			if (atoms.indexOf(at.getId()) != -1)
				return c;
		}
		return null;
	}

	public Ring getRing(Atom at) {
		for (int i = 0; i < rings.size(); i++) {
			Ring c = (Ring) rings.get(i);
			LinkedList atoms = c.getAtoms();
			if (atoms.indexOf(at.getId()) != -1)
				return c;
		}
		return null;
	}

	public Atom getAtom(int id) {
		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			if (at.getId() == id)
				return at;
		}
		return null;
	}

	public LinkedList getCloseRes(Residue res) {
		LinkedList l = new LinkedList();
		boolean found = false;

		LinkedList chains = mol.getChains();
		for (int c = 0; c < chains.size(); c++) {
			Chain ch = (Chain) chains.get(c);
			if ((ch.getId() != chain.getId()) && (!ch.getType().equals("ion"))) {
				LinkedList residues = ch.getResidues();

				for (int j = 0; j < residues.size(); j++) {
					Residue r = (Residue) residues.get(j);

					if (res.getId() != r.getId()) {
						found = false;
						LinkedList atoms1 = res.getAtoms();
						LinkedList atoms2 = r.getAtoms();
						for (int u = 0; u < atoms1.size(); u++) {
							Atom a1 = (Atom) atoms1.get(u);
							if (found)
								break;
							for (int v = 0; v < atoms2.size(); v++) {
								Atom a2 = (Atom) atoms2.get(v);
								if (isDist(a1, a2, cutoff)) {

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

	/* return the H bond list */
	public LinkedList getHBonds() {
		LinkedList l = new LinkedList();
		for (int s = 0; s < residues.size(); s++) {
			Residue res = (Residue) residues.get(s);
			for (int i = 0; i < close.size(); i++) {
				boolean bool = false;
				Residue rr = (Residue) close.get(i);
				if (res != rr) {
					for (int w = 0; w < res.atoms.size(); w++) {
						Atom a1 = (Atom) res.atoms.get(w);
						for (int j = 0; j < rr.atoms.size(); j++) {
							Atom a2 = (Atom) rr.atoms.get(j);

							if ((isHbond(a1, a2)) && (a1.distance(a2) < 3.3)) {
								boolean ok;
								if (isHBondAcceptor(a2) && isHBondDonor(a1)) {
									ok = isCompa(a1, a2);
								} else {
									ok = isCompa(a2, a1);
								}
								if (ok) {

									Hbond h;
									if (isHBondAcceptor(a2) && isHBondDonor(a1)) {
										h = new Hbond(a1, a2, rr);

									} else {
										h = new Hbond(a2, a1, rr);
									}

									l.add(h);
									bool = true;
								}
							}
						}
					}
				}

				if (bool) {
					delete(rr, copy);
				}
			}

		}

		for (int i = 0; i < l.size() - 1; i++) {
			Hbond h = (Hbond) l.get(i);
			Atom atomRes = h.getAtomRes();

			for (int j = i + 1; j < l.size(); j++) {
				Hbond h2 = (Hbond) l.get(j);
				Atom dd = h2.getDon();
				Atom aa = h2.getAcc();
				if ((atomRes.getId() == dd.getId())
						|| (atomRes.getId() == aa.getId())) {
					h.addLinked(h2);
					h2.addLinked(h);
				}
			}
		}

		for (int i = 0; i < l.size(); i++) {
			Hbond h = (Hbond) l.get(i);
			if (h.getLinked().size() == 0) {
				traitHbond(h, l);
			} else {
				traitHbond2(h, l);
			}

		}

		return l;
	}

	public LinkedList getHbondsW() {
		LinkedList res = new LinkedList();

		for (int i = 0; i < atoms.size(); i++) {
			Atom at1 = (Atom) atoms.get(i);

			for (int j = 0; j < water.size(); j++) {
				Atom at2 = (Atom) water.get(j);
				double dist1 = at1.distance(at2);

				if (dist1 < 3.3) {

					if ((isHBondDonor(at1)) || (isHBondAcceptor(at1))) {
						boolean ok = false;
						if (isHBondDonor(at1)) {
							ok = isCompa(at1, at2);
						} else {
							ok = isCompa(at2, at1);
						}
						if (ok) {

							LinkedList chains = mol.getChains();
							for (int c = 0; c < chains.size(); c++) {
								Chain ch = (Chain) chains.get(c);
								if (ch.getId() != chain.getId()) {
									LinkedList residues = ch.getResidues();

									for (int f = 0; f < residues.size(); f++) {
										Residue r = (Residue) residues.get(f);

										LinkedList atoms2 = r.getAtoms();

										for (int v = 0; v < atoms2.size(); v++) {
											Atom a2 = (Atom) atoms2.get(v);
											double dist2 = at2.distance(a2);
											if (dist2 < 3.3) {

												if ((isHBondDonor(a2))
														|| isHBondAcceptor(a2)) {
													boolean ok2 = false;
													if (isHBondDonor(a2)) {
														ok2 = isCompa(a2, at2);
													} else {
														ok2 = isCompa(at2, a2);
													}
													if (ok2) {

														Hbond h = new Hbond(
																at1, a2,
																a2.getParent(),
																dist1 + dist2);

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

			}

		}

		return res;
	}

	public boolean isCompa(Atom don, Atom acc) {

		LinkedList neig = acc.getNeig();

		for (int i = 0; i < neig.size(); i++) {
			Atom aa = (Atom) neig.get(i);

			if (getAngle3D(don, acc, aa) < 90)
				return false;
		}

		return true;
	}

	public double getAngle3D(Atom at1, Atom at2, Atom at3) {
		double angle = 0.0;
		double[] c1 = new double[3];
		double[] c2 = new double[3];
		double[] c3 = new double[3];
		c1[0] = at1.getX();
		c1[1] = at1.getY();
		c1[2] = at1.getZ();
		c2[0] = at2.getX();
		c2[1] = at2.getY();
		c2[2] = at2.getZ();
		c3[0] = at3.getX();
		c3[1] = at3.getY();
		c3[2] = at3.getZ();
		angle = angle3D(c1, c2, c3);

		return angle;
	}

	public double angle3D(double[] c1, double[] c2, double[] c3) {

		int i;
		double[] v1 = new double[3];
		double[] v2 = new double[3];
		double dot, mod1, mod2, max;
		double angle;

		for (i = 0; i < 3; i++) {
			v1[i] = c2[i] - c1[i];
			v2[i] = c3[i] - c2[i];
		}

		mod1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
		mod2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];

		if (mod1 > 0.0)
			mod1 = Math.sqrt(mod1);
		else
			return (0.0F);
		if (mod2 > 0.0)
			mod2 = Math.sqrt(mod2);
		else
			return (0.0F);

		dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

		if ((max = mod1 * mod2) < dot)
			return (0.0F);

		angle = (Math.acos(dot / max));

		angle = 180.0F - (angle * 57.29577951F);

		if (angle > 180.0F) {
			angle = angle - 180.0F;
		} else if (angle < -180.0F) {
			angle = -180.0F;
		} else if (angle < 0.0F) {
			angle = 0.0F - angle;
		}

		return angle;

	}

	public boolean isOk(Hbond h, LinkedList l) {
		Atom a1 = h.getDon();
		Atom a2 = h.getAcc();

		for (int i = 0; i < l.size(); i++) {
			Hbond h2 = (Hbond) l.get(i);
			Atom b1 = h2.getDon();
			Atom b2 = h2.getAcc();

			if ((a1.getId() == b1.getId()) || (a1.getId() == b2.getId())
					|| (a2.getId() == b1.getId()) || (a2.getId() == b2.getId()))
				return false;

		}
		return true;
	}

	public void traitHbond3(Hbond h, LinkedList l) {

		if (h.getFlag())
			return;

		LinkedList linked = h.getLinked();
		Atom aRes = h.getAtomRes();

		LinkedList ll = new LinkedList();
		ll.add(h.getAtomLig());

		for (int i = 0; i < linked.size(); i++) {
			Hbond hh = (Hbond) linked.get(i);
			hh.setFlag(true);
			Atom a = hh.getAtomLig();
			ll.add(a);
		}

		LinkedList res = barycenter(ll);
		double bx = (Double) res.get(0);
		double by = (Double) res.get(1);

		aRes.setFx(bx);
		aRes.setFy(by);

		double dist = 0;

		for (int i = 0; i < ll.size(); i++) {
			Atom aaa = (Atom) ll.get(i);

			dist = dist + aRes.distance(aaa);

		}

		dist = dist / ll.size();

		Atom temp = (Atom) ll.get(0);

		LinkedList lll = temp.getNeig();

		Atom ref = (Atom) lll.get(0);

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = aRes.getFx();
		double y1 = aRes.getFy();
		double x2 = ref.getFx();
		double y2 = ref.getFy();
		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = temp.getFx();
		double yo = temp.getFy();
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

		aRes.setFx(solx);
		aRes.setFy(soly);

	}

	public void traitHbond2(Hbond h, LinkedList l) {

		if (h.getFlag())
			return;

		LinkedList linked = h.getLinked();
		Atom aRes = h.getAtomRes();

		LinkedList ll = new LinkedList();
		ll.add(h.getAtomLig());

		for (int i = 0; i < linked.size(); i++) {
			Hbond hh = (Hbond) linked.get(i);
			hh.setFlag(true);
			Atom a = hh.getAtomLig();
			ll.add(a);
		}

		Atom at1 = (Atom) ll.get(0);
		Atom at2 = (Atom) ll.get(1);

		double dist1 = aRes.distance(at1);
		double dist2 = aRes.distance(at2);

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		double xo1 = at1.getFx();
		double yo1 = at1.getFy();
		double xo2 = at2.getFx();
		double yo2 = at2.getFy();

		double r1 = dist1 * 1.5;
		double r2 = dist2 * 1.5;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo2;
		double b2 = -2.0 * yo2;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double aa = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
		double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
		double C12 = C1 - C2;
		double bb, cc;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);

		if (bpga) {
			bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = bb * bb - 4.0 * aa * cc;

		solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
		solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta))) / aa;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			solx1 = -(b12 * solx1 + C12) / a12;
			soly2 = solx2;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		aRes.setFx(solx1);
		aRes.setFy(soly1);

		int n1 = getNumberInter(h, l);

		aRes.setFx(solx2);
		aRes.setFy(soly2);

		int n2 = getNumberInter(h, l);

		if (n1 <= n2) {
			aRes.setFx(solx1);
			aRes.setFy(soly1);
		} else {
			aRes.setFx(solx2);
			aRes.setFy(soly2);
		}

	}

	public void traitHbond(Hbond h, LinkedList l) {

		Atom aRes = h.getAtomRes();
		Atom aLig = h.getAtomLig();

		double dist = 0;

		dist = aRes.distance(aLig);

		LinkedList tmp = aLig.getNeig();
		Atom ref = (Atom) tmp.get(0);

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = aLig.getFx();
		double y1 = aLig.getFy();
		double x2 = aRes.getX();
		double y2 = aRes.getY();
		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;
		double r = dist * 1.5;

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

			aRes.setFx(sol2x);
			aRes.setFy(sol2y);
		} else {

			aRes.setFx(sol1x);
			aRes.setFy(sol1y);

		}

		double[] rot = new double[72];

		for (int s = 0; s < rot.length; s++) {
			rot[s] = 5;
		}

		double dist2 = 3.0;

		int bestInter = getNumberInter(h, l);
		double bestAtom = getNumberAtom(h, l, dist2);

		int best = -1;
		double optX = aRes.getFx();
		double optY = aRes.getFy();

		for (int v = 0; v < rot.length - 1; v++) {

			h.rotation(rot[v]);

			int bi = getNumberInter(h, l);

			if (bi < bestInter) {

				double ba = getNumberAtom(h, l, dist2);
				bestInter = bi;
				bestAtom = ba;
				best = v;
				optX = aRes.getFx();
				optY = aRes.getFy();

			} else {

				if (bi == bestInter) {

					double ba = getNumberAtom(h, l, dist2);

					if (ba > bestAtom) {
						bestInter = bi;
						bestAtom = ba;
						best = v;
						optX = aRes.getFx();
						optY = aRes.getFy();
					}
				}
			}

		}

		aRes.setFx(optX);
		aRes.setFy(optY);

		h.setFlag(true);

	}

	public boolean delete(Residue res, LinkedList l) {
		for (int i = 0; i < l.size(); i++) {
			Residue r = (Residue) l.get(i);
			if ((r.getId().equals(res.getId()))
					&& (r.getParent().getName() == res.getParent().getName())) {
				l.remove(i);
				return true;
			}
		}
		return false;
	}

	public boolean isHbond(Atom a1, Atom a2) {

		if (isHBondAcceptor(a1) && isHBondDonor(a2)) {
			return true;
		}
		if (isHBondDonor(a1) && isHBondAcceptor(a2)) {
			return true;
		}
		return false;
	}

	public final boolean isHBondAcceptor(Atom atom) {
		String type = atom.getParent().getType();

		if (type.equals("aa") || type.equals("na"))
			return atom.isHBondAcceptor();
		else {
			boolean b = atom.getParent().getMapping();
			if (b)
				return atom.isHBondAcceptor2();
			else
				return guessHbond(atom);
		}
	}

	public final boolean guessHbond(Atom atom) {

		String element = atom.getElement();
		if (element.equals("O") || element.equals("N") || element.equals("S"))
			return true;
		else
			return false;
	}

	public final boolean isHBondDonor(Atom atom) {
		String type = atom.getParent().getType();
		if (type.equals("aa") || type.equals("na"))
			return atom.isHBondDonor();
		else {
			if (atom.getParent().mapping)
				return atom.isHBondDonor2();
			else
				return guessHbond(atom);
		}

	}

	public LinkedList getCloseAtom(Residue res) {
		LinkedList l = new LinkedList();

		LinkedList chains = mol.getChains();
		for (int c = 0; c < chains.size(); c++) {
			Chain ch = (Chain) chains.get(c);
			if (ch.getId() != chain.getId()) {
				LinkedList residues = ch.getResidues();

				for (int j = 0; j < residues.size(); j++) {
					Residue r = (Residue) residues.get(j);

					if (res.getId() != r.getId()) {

						LinkedList atoms1 = res.getAtoms();
						LinkedList atoms2 = r.getAtoms();
						for (int u = 0; u < atoms1.size(); u++) {
							Atom a1 = (Atom) atoms1.get(u);

							for (int v = 0; v < atoms2.size(); v++) {
								Atom a2 = (Atom) atoms2.get(v);
								if (isDist(a1, a2, cutoff)) {
									if (!isInside(a2, l)) {
										l.add(a2);
									}
								}

							}

						}
					}
				}
			}
		}
		return l;
	}

	public boolean isInside(Atom at, LinkedList l) {
		for (int i = 0; i < l.size(); i++) {
			Atom ato = (Atom) l.get(i);
			if (at.getId() == ato.getId())
				return true;
		}
		return false;
	}

	public boolean bondExist(Bond b) {
		for (int i = 0; i < bonds.size(); i++) {
			Bond bond = (Bond) bonds.get(i);
			if (b.equals(bond))
				return true;
		}
		return false;
	}

	public boolean bondExistEx(Bond b) {
		for (int i = 0; i < explicit.size(); i++) {
			Bond bond = (Bond) explicit.get(i);
			if (b.equals(bond))
				return true;
		}
		return false;
	}

	public boolean isInside(Residue res, LinkedList l) {
		for (int i = 0; i < l.size(); i++) {
			Residue r = (Residue) l.get(i);
			if (r == res)
				return true;
		}
		return false;
	}

	public boolean isDist(Atom a1, Atom a2, double d) {
		if (a1.distance(a2) < d)
			return true;
		else
			return false;

	}

	public void updateBond(int a1, int a2, int order) {

		Bond b1 = new Bond(getAtom(a1), getAtom(a2), order);
		for (int i = 0; i < bonds.size(); i++) {
			Bond b2 = (Bond) bonds.get(i);
			if (b1.equals(b2)) {
				b2.setOrder(order);
				break;
			}
		}
	}

	public void updateBonds() {

		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);

			Graph g1 = new Graph(res);

			Graph g3 = new Graph(res.getName());
			Graph g2 = g3.withoutH();

			if (g2.isUnique()) {

				if (g2.isCompa(g1.atoms)) {

					if (g2.length >= g1.length) {

						for (int u = 0; u < g1.length; u++) {
							updateH(g1.ids[u], g3.getH(g3.getPos(g1.atoms[u])));
						}

						for (int s = 1; s < g1.length - 1; s++) {
							for (int t = s + 1; t < g1.length; t++) {

								if (g1.adja[s][t] != g2.adja[g2
										.getPos(g1.atoms[s])][g2
										.getPos(g1.atoms[t])])

								{

									updateBond(g1.ids[s], g1.ids[t], g2.adja[g2
											.getPos(g1.atoms[s])][g2
											.getPos(g1.atoms[t])]);

								}

							}
						}

						res.setMapping(true);

					} else {

					}

				} else {

				}

			} else {

			}

		}
	}

	public void updateH(int id, int h) {

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			if (at.getId() == id) {
				at.setHnumber(h);
				break;
			}
		}
	}

	public int getPos(int[] tab, int elt) {
		for (int i = 0; i < tab.length; i++) {
			if (tab[i] == elt)
				return i;
		}
		return -1;
	}

	public void print() {
		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			at.print();
		}
		for (int i = 0; i < bonds.size(); i++) {
			Bond b = (Bond) bonds.get(i);
			b.print();
		}
		System.out.println("chaines " + chaines);
	}

	public void updateAromatic() {

		for (int i = 0; i < rings.size(); i++) {
			Ring r = (Ring) rings.get(i);
			if (isAromatic(r)) {
				r.setAromatic(true);
			}
		}
		for (int i = 0; i < rings.size(); i++) {
			Ring r = (Ring) rings.get(i);
			if (isAromatic(r)) {
				r.setAromatic(true);
			}
		}
	}

	public boolean isAromatic(Ring r) {

		int simple = 0;
		int doub = 0;
		LinkedList l = r.getAtoms();
		if (l.size() != 6)
			return false;

		for (int i = 0; i < l.size() - 1; i++) {
			int id1 = (Integer) l.get(i);
			Atom at = getAtom(id1);
			if (!at.getElement().equals("C"))
				return false;
			for (int j = i + 1; j < l.size(); j++) {
				int id2 = (Integer) l.get(j);
				Bond b = getBond(id1, id2);
				if (b != null) {
					if (b.getOrder() == 1)
						simple++;
					if (b.getOrder() == 2)
						doub++;

				}
			}
		}

		if ((simple == 3) && (doub == 3))
			return true;
		LinkedList ll = r.getConjugate();
		if (ll.size() != 0) {
			for (int i = 0; i < ll.size(); i++) {
				Ring ri = (Ring) ll.get(i);
				if ((ri.isAromatic()) && (doub == 2))
					return true;
			}

		}

		return false;
	}

	public Bond getBond(int id1, int id2) {
		for (int i = 0; i < bonds.size(); i++) {
			Bond b = (Bond) bonds.get(i);
			if ((b.getFirst().getId() == id1) && (b.getSecond().getId() == id2))
				return b;
			if ((b.getFirst().getId() == id2) && (b.getSecond().getId() == id1))
				return b;
		}
		return null;
	}

	public Atom getFirst() {

		Atom res = null;

		if (chaines.size() != 0) {
			int taille = 0;
			Chaine cc = null;
			for (int i = 0; i < chaines.size(); i++) {
				Chaine chaine = (Chaine) chaines.get(0);
				if (chaine.getSize() > taille)
					cc = chaine;
			}
			LinkedList l = cc.getAtoms();

			return getAtom((Integer) l.get(l.size() - 1));
		}

		double dep = -9.9e19;
		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			if (!at.getTerm()) {
				if (at.getY() > dep) {
					res = at;
					dep = at.getY();
				}
			}
		}

		return res;

	}

	public boolean atom2Linear(Atom atom) {

		if ((atom.getType() == 1) || (atom.getType() == 3)
				|| (atom.getType() == 4) || (atom.getType() == 6))
			return true;
		LinkedList t = atom.getNeig();
		for (int i = 0; i < t.size(); i++) {
			Atom a = (Atom) t.get(i);
			if ((a.getType() == 2) && (a.getNeig().size() == 2))
				return false;
			if (a.getType() == 6)
				return false;
		}
		return true;
	}

	/* tratment of atom neighbors */
	public void traitAtom(Atom atom) {
		LinkedList neig = atom.getNeig();

		if (neig.size() == 0)
			return;

		if (neig.size() == 1) {
			traitAtom1((Atom) neig.get(0), atom);
		} else {
			if (neig.size() == 2) {
				if (!atom2Linear(atom))
					traitAtom2((Atom) neig.get(0), (Atom) neig.get(1), atom);
				else
					traitAtom2Lin((Atom) neig.get(0), (Atom) neig.get(1), atom);
			}

			else {
				if (neig.size() == 3) {
					traitAtom3((Atom) neig.get(0), (Atom) neig.get(1),
							(Atom) neig.get(2), atom);
				} else {
					if (neig.size() == 4) {
						traitAtom4((Atom) neig.get(0), (Atom) neig.get(1),
								(Atom) neig.get(2), (Atom) neig.get(3), atom);
					}

				}

			}

		}
		if ((atom.getType() == 5) || (atom.getType() == 6)) {
			atom.setFlag(true);
		}

	}

	/* traitment of a ring atom neighbors */
	public void traitAtomC(Atom atom, Atom suiv) {
		LinkedList neig = atom.getNeig();

		if (neig.size() == 3) {
			traitAtom3C(atom, suiv);
		} else {
			if (neig.size() == 4) {
				traitAtom4C(atom, suiv);
			}

		}
		if ((atom.getType() == 5) || (atom.getType() == 6)) {
			atom.setFlag(true);
		}

	}

	/* atom with one neighbor */
	private void traitAtom1(Atom trait, Atom ref) {

		if (trait.getFlag())
			return;

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = ref.getFx();
		double y1 = ref.getFy();
		double x2 = trait.getFx();
		double y2 = trait.getFy();
		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;
		double r = 1.2;

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

		double solx;
		double soly;

		if (a <= 0) {
			if (sol1y < sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		} else {
			if (sol1y > sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		}

		trait.setFx(solx);
		trait.setFy(soly);
		if ((trait.getType() == 6)) {
			trait.setFlag(true);
		}
		if ((ref.getType() == 6) || (ref.getType() == 5)
				|| (ref.getType() == 2))
			ref.setFlag(true);
	}

	private void traitAtom1C(Atom trait, Atom ref) {

		if ((trait.getFlag()) || (trait.getType() != 6))
			return;

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = ref.getFx();
		double y1 = ref.getFy();
		double x2 = trait.getFx();
		double y2 = trait.getFy();
		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;
		double r = 1.2;

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

		double solx;
		double soly;

		if (a <= 0) {
			if (sol1y < sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		} else {
			if (sol1y > sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		}

		trait.setFx(solx);
		trait.setFy(soly);
		if ((trait.getType() == 5) || (trait.getType() == 6)) {
			trait.setFlag(true);
		}
		if ((ref.getType() == 6) || (ref.getType() == 5)
				|| (ref.getType() == 2))
			ref.setFlag(true);
	}

	LinkedList atomRota;

	public LinkedList getListRota(Atom at) {

		atomRota = new LinkedList();

		LinkedList open = new LinkedList();
		LinkedList closed = new LinkedList();
		closed.add(at.getId());

		LinkedList nneig = at.getNeig();

		Atom a1 = (Atom) nneig.get(0);
		Atom a2 = (Atom) nneig.get(1);

		if (a1.getFy() > a2.getFy()) {
			open.add(a2.getId());
			closed.add(a1.getId());
		} else {
			open.add(a1.getId());
			closed.add(a2.getId());
		}

		traitRota(closed, open);

		return atomRota;

	}

	public LinkedList getListRota2(Atom at) {

		LinkedList res = new LinkedList();

		atomRota = new LinkedList();

		LinkedList open = new LinkedList();
		LinkedList closed = new LinkedList();
		closed.add(at.getId());

		LinkedList nneig = at.getNeig();

		Atom a1 = null;
		Atom a2 = null;
		Atom a3 = null;

		if (nneig.size() == 3) {
			a1 = (Atom) nneig.get(0);
			a2 = (Atom) nneig.get(1);
			a3 = (Atom) nneig.get(2);
		}

		if (nneig.size() == 2) {
			a1 = (Atom) nneig.get(0);
			a2 = (Atom) nneig.get(1);
		}

		if (a1.getFy() > a2.getFy()) {
			open.add(a2.getId());
			closed.add(a1.getId());
		} else {
			open.add(a1.getId());
			closed.add(a2.getId());
		}

		traitRota(closed, open);

		res.add(atomRota);

		atomRota = new LinkedList();

		open = new LinkedList();
		closed = new LinkedList();
		closed.add(at.getId());

		if (a3 != null) {
			open.add(a3.getId());
		}
		closed.add(a1.getId());

		traitRota(closed, open);

		res.add(atomRota);

		return res;

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

	public void Rotation(Atom pivot, LinkedList l, double angle) {

		LinkedList tr = new LinkedList();

		angle = angle * Math.PI / 180;

		for (int i = 0; i < l.size(); i++) {
			Atom at = (Atom) l.get(i);

			if ((at.getType() == 4) || (at.getType() == 3)) {
				for (int j = 0; j < rings.size(); j++) {
					Ring r = (Ring) rings.get(j);
					if (!inside(r, tr)) {
						if (r.isAromatic()) {

							LinkedList a = r.getAtoms();
							if (isInsideId(at, a)) {

								double x0 = r.getCenterX();
								double y0 = r.getCenterY();

								double OoAx = x0 - pivot.getFx();
								double OoAy = y0 - pivot.getFy();

								double rresx = pivot.getFx() + Math.cos(angle)
										* OoAx - Math.sin(angle) * OoAy;
								double rresy = pivot.getFy() + Math.sin(angle)
										* OoAx + Math.cos(angle) * OoAy;

								r.setCenterX(rresx);
								r.setCenterY(rresy);

								tr.add(r);

							}
						}

					}
				}

			}

			double OAx = at.getFx() - pivot.getFx();
			double OAy = at.getFy() - pivot.getFy();

			double resx = pivot.getFx() + Math.cos(angle) * OAx
					- Math.sin(angle) * OAy;
			double resy = pivot.getFy() + Math.sin(angle) * OAx
					+ Math.cos(angle) * OAy;

			at.setFx(resx);
			at.setFy(resy);

		}

	}

	public boolean inside(Ring r, LinkedList l) {

		for (int i = 0; i < l.size(); i++) {
			Ring ring = (Ring) l.get(i);
			if (r == ring)
				return true;

		}
		return false;

	}

	public void traitRota(LinkedList closed, LinkedList open) {

		LinkedList temp = new LinkedList();

		if (open.size() == 0) {
		} else {
			do {

				int element = (Integer) open.get(0);

				if (!isInside(getAtom(element), atomRota)) {
					atomRota.add(getAtom(element));
				}

				open.remove(0);
				temp.clear();
				Atom at = getAtom(element);

				for (int i = 0; i < at.getNeig().size(); i++) {

					Atom atom = (Atom) at.getNeig().get(i);

					temp.add(atom.getId());
				}

				for (int i = 0; i < temp.size(); i++) {

					if (estDedans((Integer) temp.get(i), closed)) {
						temp.remove(i);
						i--;
					}
				}

				closed.add(element);

				traitRota(closed, temp);

			} while (open.size() != 0);

		}

	}

	/* traitment of an atom with 2 neighbors liner way */
	public void traitAtom2Lin(Atom a1, Atom a2, Atom ref) {

		if ((a1.getFlag()) && (a2.getFlag()))
			return;

		Atom fixed = null;
		Atom var = null;

		if (a1.getFlag()) {
			fixed = a1;
			var = a2;
		} else {
			if (a2.getFlag()) {
				fixed = a2;
				var = a1;
			}
		}

		if (fixed != null) {

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = fixed.getFx();
			double y2 = fixed.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double solx;
			double soly;

			double d1 = Math.sqrt((fixed.getFx() - sol1x)
					* (fixed.getFx() - sol1x) + (fixed.getFy() - sol1y)
					* (fixed.getFy() - sol1y));
			double d2 = Math.sqrt((fixed.getFx() - sol2x)
					* (fixed.getFx() - sol2x) + (fixed.getFy() - sol2y)
					* (fixed.getFy() - sol2y));

			if (d1 > d2) {
				solx = sol1x;
				soly = sol1y;
			} else {
				solx = sol2x;
				soly = sol2y;
			}

			var.setFx(solx);
			var.setFy(soly);
			if ((var.getType() == 6)) {
				var.setFlag(true);
			}
		}

		else {
			Atom ref2 = null;
			Atom var2 = null;

			if ((a1.getTerm()) && (a2.getTerm())) {
				ref2 = a1;
				var2 = a2;
			} else {
				if (a1.getTerm()) {
					var2 = a1;
					ref2 = a2;
				} else {
					if (a2.getTerm()) {
						var2 = a2;
						ref2 = a1;
					} else {
						ref2 = a1;
						var2 = a2;
					}
				}
			}

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = ref2.getFx();
			double y2 = ref2.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double d1 = Math.sqrt((ref2.getFx() - sol1x)
					* (ref2.getFx() - sol1x) + (ref2.getFy() - sol1y)
					* (ref2.getFy() - sol1y));
			double d2 = Math.sqrt((var2.getFx() - sol1x)
					* (var2.getFx() - sol1x) + (var2.getFy() - sol1y)
					* (var2.getFy() - sol1y));

			if (d1 < d2) {
				ref2.setFx(sol1x);
				ref2.setFy(sol1y);
				if ((ref2.getType() == 6)) {
					ref2.setFlag(true);
				}
				var2.setFx(sol2x);
				var2.setFy(sol2y);
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}
			} else {
				ref2.setFx(sol2x);
				ref2.setFy(sol2y);
				if ((ref2.getType() == 6)) {
					ref2.setFlag(true);
				}
				var2.setFx(sol1x);
				var2.setFy(sol1y);
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}
			}

		}
		if ((ref.getType() == 6) || (ref.getType() == 5)
				|| (ref.getType() == 2))
			ref.setFlag(true);
	}

	/* traitment of an atom with 2 neighbors zig-zag way */
	public void traitAtom2(Atom a1, Atom a2, Atom ref) {

		if ((ref.getType() == 5) || (ref.getType() == 2)) {
			if ((!a1.getTerm()) && (!a2.getTerm()))
				ref.setPivot(true);
		}

		if ((a1.getFlag()) && (a2.getFlag()))
			return;

		Atom fixed = null;
		Atom var = null;

		if (a1.getFlag()) {
			fixed = a1;
			var = a2;
		} else {
			if (a2.getFlag()) {
				fixed = a2;
				var = a1;
			}
		}

		if (fixed != null) {

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = fixed.getFx();
			double y2 = fixed.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double solx;
			double soly;

			double d1 = Math.sqrt((fixed.getFx() - sol1x)
					* (fixed.getFx() - sol1x) + (fixed.getFy() - sol1y)
					* (fixed.getFy() - sol1y));
			double d2 = Math.sqrt((fixed.getFx() - sol2x)
					* (fixed.getFx() - sol2x) + (fixed.getFy() - sol2y)
					* (fixed.getFy() - sol2y));

			if (d1 > d2) {
				solx = sol1x;
				soly = sol1y;
			} else {
				solx = sol2x;
				soly = sol2y;
			}

			double OAx = solx - ref.getFx();
			double OAy = soly - ref.getFy();

			double resx = ref.getFx() + Math.cos(1.04) * OAx - Math.sin(1.04)
					* OAy;
			double resy = ref.getFy() + Math.sin(1.04) * OAx + Math.cos(1.04)
					* OAy;

			double resx2 = ref.getFx() + Math.cos(-1.04) * OAx
					- Math.sin(-1.04) * OAy;
			double resy2 = ref.getFy() + Math.sin(-1.04) * OAx
					+ Math.cos(-1.04) * OAy;

			double conf1 = getAtomConflit(resx, resy);
			double conf2 = getAtomConflit(resx2, resy2);

			if (conf1 != conf2) {

				if (conf2 > conf1) {
					var.setFx(resx);
					var.setFy(resy);
				} else {
					var.setFx(resx2);
					var.setFy(resy2);
				}
				if ((var.getType() == 6)) {
					var.setFlag(true);
				}

			} else {

				LinkedList t = fixed.getNeig();
				Atom temp = null;
				for (int i = 0; i < t.size(); i++) {
					Atom aaa = (Atom) t.get(i);
					if (aaa.getId() != ref.getId()) {
						temp = aaa;
						break;
					}
				}
				boolean ff = false;
				if (temp != null) {
					if (angle(temp, fixed, ref) > 0) {
						if (angleNew(fixed.getFx(), fixed.getFy(), ref, resx,
								resy) < 0) {
							var.setFx(resx);
							var.setFy(resy);
							ff = true;
						} else {
							var.setFx(resx2);
							var.setFy(resy2);
							ff = true;
						}
					} else {
						if (angleNew(fixed.getFx(), fixed.getFy(), ref, resx,
								resy) > 0) {
							var.setFx(resx);
							var.setFy(resy);
							ff = true;
						} else {
							var.setFx(resx2);
							var.setFy(resy2);
							ff = true;
						}

					}
				}

				if (!ff) {
					var.setFx(resx);
					var.setFy(resy);
				}

				if ((var.getType() == 6)) {
					var.setFlag(true);
				}
			}
		}

		else {
			Atom ref2 = null;
			Atom var2 = null;

			if ((a1.getTerm()) && (a2.getTerm())) {
				ref2 = a1;
				var2 = a2;
			} else {
				if (a1.getTerm()) {
					var2 = a1;
					ref2 = a2;
				} else {
					if (a2.getTerm()) {
						var2 = a2;
						ref2 = a1;
					} else {
						ref2 = a1;
						var2 = a2;
					}
				}
			}

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = ref2.getFx();
			double y2 = ref2.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double d1 = Math.sqrt((ref2.getFx() - sol1x)
					* (ref2.getFx() - sol1x) + (ref2.getFy() - sol1y)
					* (ref2.getFy() - sol1y));
			double d2 = Math.sqrt((var2.getFx() - sol1x)
					* (var2.getFx() - sol1x) + (var2.getFy() - sol1y)
					* (var2.getFy() - sol1y));

			if (d1 < d2) {
				ref2.setFx(sol1x);
				ref2.setFy(sol1y);
				if ((ref2.getType() == 6)) {
					ref2.setFlag(true);
				}
				var2.setFx(sol2x);
				var2.setFy(sol2y);
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}
			} else {
				ref2.setFx(sol2x);
				ref2.setFy(sol2y);
				if ((ref2.getType() == 6)) {
					ref2.setFlag(true);
				}
				var2.setFx(sol1x);
				var2.setFy(sol1y);
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}
			}

			double OAx = var2.getFx() - ref.getFx();
			double OAy = var2.getFy() - ref.getFy();

			double resx = ref.getFx() + Math.cos(1.04) * OAx - Math.sin(1.04)
					* OAy;
			double resy = ref.getFy() + Math.sin(1.04) * OAx + Math.cos(1.04)
					* OAy;

			double resx2 = ref.getFx() + Math.cos(-1.04) * OAx
					- Math.sin(-1.04) * OAy;
			double resy2 = ref.getFy() + Math.sin(-1.04) * OAx
					+ Math.cos(-1.04) * OAy;

			double conf1 = getAtomConflit(resx, resy);
			double conf2 = getAtomConflit(resx2, resy2);

			if (conf1 != conf2) {

				if (conf2 > conf1) {
					var2.setFx(resx);
					var2.setFy(resy);
				} else {
					var2.setFx(resx2);
					var2.setFy(resy2);
				}
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}

			} else {

				double d3 = Math.sqrt((var2.getX() - resx)
						* (var2.getX() - resx) + (var2.getY() - resy)
						* (var2.getY() - resy));
				double d4 = Math.sqrt((var2.getX() - resx2)
						* (var2.getX() - resx2) + (var2.getY() - resy2)
						* (var2.getY() - resy2));

				if (d3 < d4) {
					var2.setFx(resx);
					var2.setFy(resy);
				} else {
					var2.setFx(resx2);
					var2.setFy(resy2);
				}
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}
			}

		}
		if ((ref.getType() == 6) || (ref.getType() == 5)
				|| (ref.getType() == 2))
			ref.setFlag(true);
	}

	/* traitment of an atom with 2 neighbors */
	public void traitAtom2Bis(Atom a1, Atom a2, Atom ref) {

		if ((ref.getType() == 5) || (ref.getType() == 2)) {
			if ((!a1.getTerm()) && (!a2.getTerm()))
				ref.setPivot(true);
		}

		if ((a1.getFlag()) && (a2.getFlag()))
			return;

		Atom fixed = null;
		Atom var = null;

		if (a1.getFlag()) {
			fixed = a1;
			var = a2;
		} else {
			if (a2.getFlag()) {
				fixed = a2;
				var = a1;
			}
		}

		if (fixed != null) {

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = fixed.getFx();
			double y2 = fixed.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double solx;
			double soly;

			double d1 = Math.sqrt((fixed.getFx() - sol1x)
					* (fixed.getFx() - sol1x) + (fixed.getFy() - sol1y)
					* (fixed.getFy() - sol1y));
			double d2 = Math.sqrt((fixed.getFx() - sol2x)
					* (fixed.getFx() - sol2x) + (fixed.getFy() - sol2y)
					* (fixed.getFy() - sol2y));

			if (d1 > d2) {
				solx = sol1x;
				soly = sol1y;
			} else {
				solx = sol2x;
				soly = sol2y;
			}

			double OAx = solx - ref.getFx();
			double OAy = soly - ref.getFy();

			double resx = ref.getFx() + Math.cos(1.04) * OAx - Math.sin(1.04)
					* OAy;
			double resy = ref.getFy() + Math.sin(1.04) * OAx + Math.cos(1.04)
					* OAy;

			double resx2 = ref.getFx() + Math.cos(-1.04) * OAx
					- Math.sin(-1.04) * OAy;
			double resy2 = ref.getFy() + Math.sin(-1.04) * OAx
					+ Math.cos(-1.04) * OAy;

			double conf1 = getAtomConflit(resx, resy);
			double conf2 = getAtomConflit(resx2, resy2);

			if (conf1 != conf2) {

				if (conf2 > conf1) {
					var.setFx(resx);
					var.setFy(resy);
				} else {
					var.setFx(resx2);
					var.setFy(resy2);
				}
				if ((var.getType() == 6)) {
					var.setFlag(true);
				}

			} else {

				double d3 = Math.sqrt((var.getX() - resx) * (var.getX() - resx)
						+ (var.getY() - resy) * (var.getY() - resy));
				double d4 = Math.sqrt((var.getX() - resx2)
						* (var.getX() - resx2) + (var.getY() - resy2)
						* (var.getY() - resy2));

				if (d3 < d4) {
					var.setFx(resx);
					var.setFy(resy);
				} else {
					var.setFx(resx2);
					var.setFy(resy2);
				}
				if ((var.getType() == 6)) {
					var.setFlag(true);
				}
			}
		}

		else {
			Atom ref2 = null;
			Atom var2 = null;

			if ((a1.getTerm()) && (a2.getTerm())) {
				ref2 = a1;
				var2 = a2;
			} else {
				if (a1.getTerm()) {
					var2 = a1;
					ref2 = a2;
				} else {
					if (a2.getTerm()) {
						var2 = a2;
						ref2 = a1;
					} else {
						ref2 = a1;
						var2 = a2;
					}
				}
			}

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = ref2.getFx();
			double y2 = ref2.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double d1 = Math.sqrt((ref2.getFx() - sol1x)
					* (ref2.getFx() - sol1x) + (ref2.getFy() - sol1y)
					* (ref2.getFy() - sol1y));
			double d2 = Math.sqrt((var2.getFx() - sol1x)
					* (var2.getFx() - sol1x) + (var2.getFy() - sol1y)
					* (var2.getFy() - sol1y));

			if (d1 < d2) {
				ref2.setFx(sol1x);
				ref2.setFy(sol1y);
				if ((ref2.getType() == 6)) {
					ref2.setFlag(true);
				}
				var2.setFx(sol2x);
				var2.setFy(sol2y);
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}
			} else {
				ref2.setFx(sol2x);
				ref2.setFy(sol2y);
				if ((ref2.getType() == 6)) {
					ref2.setFlag(true);
				}
				var2.setFx(sol1x);
				var2.setFy(sol1y);
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}
			}

			double OAx = var2.getFx() - ref.getFx();
			double OAy = var2.getFy() - ref.getFy();

			double resx = ref.getFx() + Math.cos(1.04) * OAx - Math.sin(1.04)
					* OAy;
			double resy = ref.getFy() + Math.sin(1.04) * OAx + Math.cos(1.04)
					* OAy;

			double resx2 = ref.getFx() + Math.cos(-1.04) * OAx
					- Math.sin(-1.04) * OAy;
			double resy2 = ref.getFy() + Math.sin(-1.04) * OAx
					+ Math.cos(-1.04) * OAy;

			double conf1 = getAtomConflit(resx, resy);
			double conf2 = getAtomConflit(resx2, resy2);

			if (conf1 != conf2) {

				if (conf2 > conf1) {
					var2.setFx(resx);
					var2.setFy(resy);
				} else {
					var2.setFx(resx2);
					var2.setFy(resy2);
				}
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}

			} else {

				double d3 = Math.sqrt((var2.getX() - resx)
						* (var2.getX() - resx) + (var2.getY() - resy)
						* (var2.getY() - resy));
				double d4 = Math.sqrt((var2.getX() - resx2)
						* (var2.getX() - resx2) + (var2.getY() - resy2)
						* (var2.getY() - resy2));

				if (d3 < d4) {
					var2.setFx(resx);
					var2.setFy(resy);
				} else {
					var2.setFx(resx2);
					var2.setFy(resy2);
				}
				if ((var2.getType() == 6)) {
					var2.setFlag(true);
				}
			}

		}
		if ((ref.getType() == 6) || (ref.getType() == 5)
				|| (ref.getType() == 2))
			ref.setFlag(true);
	}

	public double getAtomConflit(double x, double y) {

		double res = 0;

		for (int i = 0; i < atoms.size(); i++) {
			Atom a = (Atom) atoms.get(i);
			if (a.getFlag()) {
				double dis = Math.sqrt((x - a.getFx()) * (x - a.getFx())
						+ (y - a.getFy()) * (x - a.getFy()));
				if (dis < 4) {
					res = res + dis;
				}

			}

		}
		return res;

	}

	/* traitment of an atom with 3 neighbors */
	public void traitAtom3(Atom a1, Atom a2, Atom a3, Atom ref) {

		if ((ref.getType() == 5) || (ref.getType() == 2)) {
			if (!a1.getTerm())
				ref.setPivot2(true);
		}

		if ((a1.getFlag()) && (a2.getFlag()) && (a3.getFlag()))
			return;
		int ct = 0;
		if (a1.getFlag())
			ct = ct + 1;
		if (a2.getFlag())
			ct = ct + 1;
		if (a3.getFlag())
			ct = ct + 1;

		if (ct == 0) {

			if (!a1.getTerm()) {
				twoCercles(a1, ref, a2, a3);
			} else {
				if (!a2.getTerm()) {
					twoCercles(a2, ref, a1, a3);
				} else {
					if (!a3.getTerm()) {
						twoCercles(a3, ref, a1, a2);
					} else {
						twoCercles(a1, ref, a2, a3);
					}

				}
			}

		}

		else {
			if ((ct == 1) || (ct == 2)) {

				if (a1.getFlag()) {
					twoCercles(a1, ref, a2, a3);
				} else {
					if (a2.getFlag()) {
						twoCercles(a2, ref, a1, a3);
					} else {
						twoCercles(a3, ref, a1, a2);
					}
				}
			}

		}
		if ((ref.getType() == 6) || (ref.getType() == 5)
				|| (ref.getType() == 2))
			ref.setFlag(true);
	}

	/* traitment of an atom with 4 neighbors */
	public void traitAtom4(Atom a1, Atom a2, Atom a3, Atom a4, Atom ref) {

		Atom[] tab = new Atom[4];

		tab[0] = a1;
		tab[1] = a2;
		tab[2] = a3;
		tab[3] = a4;

		Atom[] newTab = rigthOrder(ref, tab);

		int ct = 0;
		if (a1.getFlag()) {
			ct = ct + 1;
		}
		if (a2.getFlag()) {
			ct = ct + 1;
		}
		if (a3.getFlag()) {
			ct = ct + 1;
		}
		if (a4.getFlag()) {
			ct = ct + 1;
		}

		if (ct == 4)
			return;
		else {
			if (ct == 0) {
				twoCercles(ref, newTab);
			} else {
				twoCercles2(ref, newTab);
			}
		}
		if ((ref.getType() == 6) || (ref.getType() == 5)
				|| (ref.getType() == 2))
			ref.setFlag(true);
	}

	/* traitment of a ring atom with 3 neighbors */
	public void traitAtom3C(Atom ref, Atom suiv) {

		LinkedList neig = ref.getNeig();
		Atom a1 = null;
		Atom a2 = null;

		Atom att1 = (Atom) neig.get(0);
		Atom att2 = (Atom) neig.get(1);
		Atom att3 = (Atom) neig.get(2);

		if (att1.getId() != suiv.getId()) {
			a1 = att1;
			if (att2.getId() != suiv.getId()) {
				a2 = att2;
			} else {
				a2 = att3;
			}
		} else {
			a1 = att2;
			a2 = att3;
		}

		if ((a1.getFlag()) && (a2.getFlag()))
			return;

		twoCercles2(ref, suiv, a1, a2);

		if ((ref.getType() == 6) || (ref.getType() == 5)
				|| (ref.getType() == 2))
			ref.setFlag(true);
	}

	/* traitment of a ring atom with 4 neighbors */
	public void traitAtom4C(Atom ref, Atom suiv) {

		LinkedList neig = ref.getNeig();

		Atom a1 = null;
		Atom a2 = null;
		Atom a3 = null;

		Atom att1 = (Atom) neig.get(0);
		Atom att2 = (Atom) neig.get(1);
		Atom att3 = (Atom) neig.get(2);
		Atom att4 = (Atom) neig.get(3);

		if (att1.getId() != suiv.getId()) {
			a1 = att1;
			if (att2.getId() != suiv.getId()) {
				a2 = att2;
				if (att3.getId() != suiv.getId())
					a3 = att3;
				else
					a3 = att4;
			} else {
				a2 = att3;
				a3 = att4;
			}

		} else {
			a1 = att2;
			a2 = att3;
			a3 = att4;
		}

		twoCercles2(ref, suiv, a1, a2, a3);

		if ((ref.getType() == 6) || (ref.getType() == 5)
				|| (ref.getType() == 2))
			ref.setFlag(true);
	}

	public void twoCercles(Atom fixed, Atom ref, Atom at1, Atom at2) {

		Atom temp = fixed.clone();

		if ((!fixed.getFlag()) || (ref.distanceF(fixed) - 1.2 > 0.1)) {

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = fixed.getFx();
			double y2 = fixed.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double d1 = Math.sqrt((fixed.getFx() - sol1x)
					* (fixed.getFx() - sol1x) + (fixed.getFy() - sol1y)
					* (fixed.getFy() - sol1y));
			double d2 = Math.sqrt((fixed.getFx() - sol2x)
					* (fixed.getFx() - sol2x) + (fixed.getFy() - sol2y)
					* (fixed.getFy() - sol2y));

			if (d1 < d2) {
				fixed.setFx(sol1x);
				fixed.setFy(sol1y);
				if ((fixed.getType() == 5) || (fixed.getType() == 6)) {
					fixed.setFlag(true);
				}
			} else {
				fixed.setFx(sol2x);
				fixed.setFy(sol2y);
				if ((fixed.getType() == 5) || (fixed.getType() == 6)) {
					fixed.setFlag(true);
				}
			}

		}

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		double xo1 = ref.getFx();
		double yo1 = ref.getFy();
		double xo2 = fixed.getFx();
		double yo2 = fixed.getFy();

		double r1 = 1.2;
		double r2 = 2.08;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo2;
		double b2 = -2.0 * yo2;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double a = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
		double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
		double C12 = C1 - C2;
		double b, c;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);

		if (bpga) {
			b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = b * b - 4.0 * a * c;

		solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
		solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			solx1 = -(b12 * solx1 + C12) / a12;
			soly2 = solx2;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		if (Math.abs(at1.getFx() - at2.getFx()) > Math.abs(at1.getFy()
				- at2.getFy())) {
			if (solx1 > solx2) {
				if (!at1.getFlag()) {
					at1.setFx(solx1);
					at1.setFy(soly1);
					if ((at1.getType() == 5) || (at1.getType() == 6)) {
						at1.setFlag(true);
					}
				}
				if (!at2.getFlag()) {
					at2.setFx(solx2);
					at2.setFy(soly2);
					if ((at2.getType() == 5) || (at2.getType() == 6)) {
						at2.setFlag(true);
					}
				}
			} else {
				if (!at1.getFlag()) {
					at1.setFx(solx2);
					at1.setFy(soly2);
					if ((at1.getType() == 5) || (at1.getType() == 6)) {
						at1.setFlag(true);
					}
				}
				if (!at2.getFlag()) {
					at2.setFx(solx1);
					at2.setFy(soly1);
					if ((at2.getType() == 5) || (at2.getType() == 6)) {
						at2.setFlag(true);
					}
				}
			}
		} else {

			if (soly1 > soly2) {
				if (!at1.getFlag()) {
					at1.setFx(solx1);
					at1.setFy(soly1);
					if ((at1.getType() == 5) || (at1.getType() == 6)) {
						at1.setFlag(true);
					}
				}
				if (!at2.getFlag()) {
					at2.setFx(solx2);
					at2.setFy(soly2);
					if ((at2.getType() == 5) || (at2.getType() == 6)) {
						at2.setFlag(true);
					}
				}
			} else {
				if (!at1.getFlag()) {
					at1.setFx(solx2);
					at1.setFy(soly2);
					if ((at1.getType() == 5) || (at1.getType() == 6)) {
						at1.setFlag(true);
					}
				}
				if (!at2.getFlag()) {
					at2.setFx(solx1);
					at2.setFy(soly1);
					if ((at2.getType() == 5) || (at2.getType() == 6)) {
						at2.setFlag(true);
					}
				}
			}

		}

	}

	public void twoCercles(Atom ref, Atom[] tab) {

		Atom top = tab[0];
		Atom rigth = tab[1];
		Atom bottom = tab[2];
		Atom left = tab[3];

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = ref.getFx();
		double y1 = ref.getFy();
		double x2 = top.getFx();
		double y2 = top.getFy();
		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;
		double r = 1.2;

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

		double d1 = Math.sqrt((top.getFx() - sol1x) * (top.getFx() - sol1x)
				+ (top.getFy() - sol1y) * (top.getFy() - sol1y));
		double d2 = Math.sqrt((top.getFx() - sol2x) * (top.getFx() - sol2x)
				+ (top.getFy() - sol2y) * (top.getFy() - sol2y));

		if (d1 < d2) {
			if (!top.getFlag()) {
				top.setFx(sol1x);
				top.setFy(sol1y);
				if ((top.getType() == 6)) {
					top.setFlag(true);
				}
			}
			if (!bottom.getFlag()) {
				bottom.setFx(sol2x);
				bottom.setFy(sol2y);
				if ((bottom.getType() == 6)) {
					bottom.setFlag(true);
				}
			}
		} else {
			if (!top.getFlag()) {
				top.setFx(sol2x);
				top.setFy(sol2y);
				if ((top.getType() == 6)) {
					top.setFlag(true);
				}
			}
			if (!bottom.getFlag()) {
				bottom.setFx(sol1x);
				bottom.setFy(sol1y);
				if ((bottom.getType() == 6)) {
					bottom.setFlag(true);
				}
			}
		}

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		double xo1 = ref.getFx();
		double yo1 = ref.getFy();
		double xo2 = top.getFx();
		double yo2 = top.getFy();

		double r1 = 1.2;
		double r2 = 1.7;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo2;
		double b2 = -2.0 * yo2;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double aa = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
		double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
		double C12 = C1 - C2;
		double bb, cc;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);

		if (bpga) {
			bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = bb * bb - 4.0 * aa * cc;

		solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
		solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta))) / aa;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			solx1 = -(b12 * solx1 + C12) / a12;
			soly2 = solx2;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		if (solx1 > solx2) {
			if (!rigth.getFlag()) {
				rigth.setFx(solx1);
				rigth.setFy(soly1);
				if ((rigth.getType() == 6)) {
					rigth.setFlag(true);
				}
			}
			if (!left.getFlag()) {
				left.setFx(solx2);
				left.setFy(soly2);
				if ((left.getType() == 6)) {
					left.setFlag(true);
				}
			}
		} else {
			if (!rigth.getFlag()) {
				rigth.setFx(solx2);
				rigth.setFy(soly2);
				if ((rigth.getType() == 6)) {
					rigth.setFlag(true);
				}
			}
			if (!left.getFlag()) {
				left.setFx(solx1);
				left.setFy(soly1);
				if ((left.getType() == 6)) {
					left.setFlag(true);
				}
			}

		}

	}

	public Atom[] rigthOrder(Atom ref, Atom[] tab) {

		Atom[] res = new Atom[4];

		double sx1, sx2, sx3, sx4;
		double sy1, sy2, sy3, sy4;

		double xo = ref.getFx();
		double yo = ref.getFy();

		double min = 200;

		for (int i = 0; i < tab.length; i++) {
			double ttt = Math.sqrt((xo - tab[i].getFx())
					* (xo - tab[i].getFx()) + (yo - tab[i].getFy())
					* (yo - tab[i].getFy()));

			if (ttt < min) {
				min = ttt;
			}
		}

		double r = min;

		LinkedList l1 = interDC(xo, yo, r, tab[0].getFx(), tab[0].getFy());
		sx1 = (Double) l1.get(0);
		sy1 = (Double) l1.get(1);
		LinkedList l2 = interDC(xo, yo, r, tab[1].getFx(), tab[1].getFy());
		sx2 = (Double) l2.get(0);
		sy2 = (Double) l2.get(1);
		LinkedList l3 = interDC(xo, yo, r, tab[2].getFx(), tab[2].getFy());
		sx3 = (Double) l3.get(0);
		sy3 = (Double) l3.get(1);
		LinkedList l4 = interDC(xo, yo, r, tab[3].getFx(), tab[3].getFy());
		sx4 = (Double) l4.get(0);
		sy4 = (Double) l4.get(1);

		double[] xs = new double[4];

		xs[0] = sx1;
		xs[1] = sx2;
		xs[2] = sx3;
		xs[3] = sx4;

		double[] ys = new double[4];

		ys[0] = sy1;
		ys[1] = sy2;
		ys[2] = sy3;
		ys[3] = sy4;

		LinkedList inter = new LinkedList();

		Atom top = null;
		Atom bottom = null;
		Atom left = null;
		Atom rigth = null;

		double ymax = -200;
		int index = -1;
		for (int i = 0; i < ys.length; i++) {
			if (ys[i] > ymax) {
				if (inter.indexOf(i) == -1) {
					ymax = ys[i];
					index = i;
				}
			}
		}
		top = tab[index];
		inter.add(index);

		double ymin = 200;
		int index2 = -1;
		for (int i = 0; i < ys.length; i++) {
			if (ys[i] < ymin) {
				if (inter.indexOf(i) == -1) {
					ymin = ys[i];
					index2 = i;
				}
			}
		}
		bottom = tab[index2];
		inter.add(index2);

		double xmax = -200;
		int index3 = -1;
		for (int i = 0; i < xs.length; i++) {
			if (xs[i] > xmax) {
				if (inter.indexOf(i) == -1) {
					xmax = xs[i];
					index3 = i;
				}
			}
		}
		rigth = tab[index3];
		inter.add(index3);

		double xmin = 200;
		int index4 = -1;
		for (int i = 0; i < xs.length; i++) {
			if (xs[i] < xmin) {
				if (inter.indexOf(i) == -1) {
					xmin = ys[i];
					index4 = i;
				}
			}
		}
		left = tab[index4];
		inter.add(index4);

		res[0] = top;
		res[1] = rigth;
		res[2] = bottom;
		res[3] = left;

		return res;

	}

	public void twoCercles2(Atom ref, Atom[] tab) {

		Atom top = tab[0];
		Atom rigth = tab[1];
		Atom bottom = tab[2];
		Atom left = tab[3];

		if (top.getFlag()) {
			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = top.getFx();
			double y2 = top.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double d1 = Math.sqrt((top.getFx() - sol1x) * (top.getFx() - sol1x)
					+ (top.getFy() - sol1y) * (top.getFy() - sol1y));
			double d2 = Math.sqrt((top.getFx() - sol2x) * (top.getFx() - sol2x)
					+ (top.getFy() - sol2y) * (top.getFy() - sol2y));

			if (d1 < d2) {
				if (!top.getFlag()) {
					top.setFx(sol1x);
					top.setFy(sol1y);
					if ((top.getType() == 6)) {
						top.setFlag(true);
					}
				}
				if (!bottom.getFlag()) {
					bottom.setFx(sol2x);
					bottom.setFy(sol2y);
					if ((bottom.getType() == 6)) {
						bottom.setFlag(true);
					}
				}
			} else {
				if (!top.getFlag()) {
					top.setFx(sol2x);
					top.setFy(sol2y);
					if ((top.getType() == 6)) {
						top.setFlag(true);
					}
				}
				if (!bottom.getFlag()) {
					bottom.setFx(sol1x);
					bottom.setFy(sol1y);
					if ((bottom.getType() == 6)) {
						bottom.setFlag(true);
					}
				}
			}

			double solx1 = 0;
			double soly1 = 0;
			double solx2 = 0;
			double soly2 = 0;

			double xo1 = ref.getFx();
			double yo1 = ref.getFy();
			double xo2 = top.getFx();
			double yo2 = top.getFy();

			double r1 = 1.2;
			double r2 = 1.7;

			double a1 = -2.0 * xo1;
			double b1 = -2.0 * yo1;
			double a2 = -2.0 * xo2;
			double b2 = -2.0 * yo2;
			double a12 = a1 - a2;
			double b12 = b1 - b2;
			double aa = a12 * a12 + b12 * b12;

			double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
			double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
			double C12 = C1 - C2;
			double bb, cc;
			boolean bpga = Math.abs(b12) >= Math.abs(a12);

			if (bpga) {
				bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
				cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
			} else {
				bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
				cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
			}
			double delta = bb * bb - 4.0 * aa * cc;

			solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
			solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta))) / aa;

			if (bpga) {
				soly1 = -(a12 * solx1 + C12) / b12;
				soly2 = -(a12 * solx2 + C12) / b12;
			} else {
				soly1 = solx1;
				solx1 = -(b12 * solx1 + C12) / a12;
				soly2 = solx2;
				solx2 = -(b12 * solx2 + C12) / a12;
			}

			if (solx1 > solx2) {
				if (!rigth.getFlag()) {
					rigth.setFx(solx1);
					rigth.setFy(soly1);
					if ((rigth.getType() == 6)) {
						rigth.setFlag(true);
					}
				}
				if (!left.getFlag()) {
					left.setFx(solx2);
					left.setFy(soly2);
					if ((left.getType() == 6)) {
						left.setFlag(true);
					}
				}
			} else {
				if (!rigth.getFlag()) {
					rigth.setFx(solx2);
					rigth.setFy(soly2);
					if ((rigth.getType() == 6)) {
						rigth.setFlag(true);
					}
				}
				if (!left.getFlag()) {
					left.setFx(solx1);
					left.setFy(soly1);
					if ((left.getType() == 6)) {
						left.setFlag(true);
					}
				}

			}

		} else {
			if (bottom.getFlag()) {

				double sol1x = 0;
				double sol1y = 0;
				double sol2x = 0;
				double sol2y = 0;

				double x1 = ref.getFx();
				double y1 = ref.getFy();
				double x2 = bottom.getFx();
				double y2 = bottom.getFy();
				double a = y2 - y1;
				double b = x1 - x2;
				double c = y1 * x2 - x1 * y2;

				double xo = x1;
				double yo = y1;
				double r = 1.2;

				if (b != 0.0) {
					double u = a * (c + b * yo) - b * b * xo;
					double v = a * a + b * b;
					double w = c + b * yo;
					double deltap = u * u - v
							* (b * b * (xo * xo - r * r) + w * w);
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

				double d1 = Math.sqrt((bottom.getFx() - sol1x)
						* (bottom.getFx() - sol1x) + (bottom.getFy() - sol1y)
						* (bottom.getFy() - sol1y));
				double d2 = Math.sqrt((bottom.getFx() - sol2x)
						* (bottom.getFx() - sol2x) + (bottom.getFy() - sol2y)
						* (bottom.getFy() - sol2y));

				if (d1 < d2) {

					if (!top.getFlag()) {
						top.setFx(sol2x);
						top.setFy(sol2y);
						if ((top.getType() == 6)) {
							top.setFlag(true);
						}
					}
				} else {

					if (!top.getFlag()) {
						top.setFx(sol1x);
						top.setFy(sol1y);
						if ((top.getType() == 6)) {
							top.setFlag(true);
						}
					}
				}

				double solx1 = 0;
				double soly1 = 0;
				double solx2 = 0;
				double soly2 = 0;

				double xo1 = ref.getFx();
				double yo1 = ref.getFy();
				double xo2 = bottom.getFx();
				double yo2 = bottom.getFy();

				double r1 = 1.2;
				double r2 = 1.7;

				double a1 = -2.0 * xo1;
				double b1 = -2.0 * yo1;
				double a2 = -2.0 * xo2;
				double b2 = -2.0 * yo2;
				double a12 = a1 - a2;
				double b12 = b1 - b2;
				double aa = a12 * a12 + b12 * b12;

				double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
				double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
				double C12 = C1 - C2;
				double bb, cc;
				boolean bpga = Math.abs(b12) >= Math.abs(a12);

				if (bpga) {
					bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
					cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
				} else {
					bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
					cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
				}
				double delta = bb * bb - 4.0 * aa * cc;

				solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
				solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta))) / aa;

				if (bpga) {
					soly1 = -(a12 * solx1 + C12) / b12;
					soly2 = -(a12 * solx2 + C12) / b12;
				} else {
					soly1 = solx1;
					solx1 = -(b12 * solx1 + C12) / a12;
					soly2 = solx2;
					solx2 = -(b12 * solx2 + C12) / a12;
				}

				if (solx1 > solx2) {
					if (!rigth.getFlag()) {
						rigth.setFx(solx1);
						rigth.setFy(soly1);
						if ((rigth.getType() == 6)) {
							rigth.setFlag(true);
						}
					}
					if (!left.getFlag()) {
						left.setFx(solx2);
						left.setFy(soly2);
						if ((left.getType() == 6)) {
							left.setFlag(true);
						}
					}
				} else {
					if (!rigth.getFlag()) {
						rigth.setFx(solx2);
						rigth.setFy(soly2);
						if ((rigth.getType() == 6)) {
							rigth.setFlag(true);
						}
					}
					if (!left.getFlag()) {
						left.setFx(solx1);
						left.setFy(soly1);
						if ((left.getType() == 6)) {
							left.setFlag(true);
						}
					}

				}

			} else {
				if (rigth.getFlag()) {

					double sol1x = 0;
					double sol1y = 0;
					double sol2x = 0;
					double sol2y = 0;

					double x1 = ref.getFx();
					double y1 = ref.getFy();
					double x2 = rigth.getFx();
					double y2 = rigth.getFy();
					double a = y2 - y1;
					double b = x1 - x2;
					double c = y1 * x2 - x1 * y2;

					double xo = x1;
					double yo = y1;
					double r = 1.2;

					if (b != 0.0) {
						double u = a * (c + b * yo) - b * b * xo;
						double v = a * a + b * b;
						double w = c + b * yo;
						double deltap = u * u - v
								* (b * b * (xo * xo - r * r) + w * w);
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

					double d1 = Math.sqrt((rigth.getFx() - sol1x)
							* (rigth.getFx() - sol1x) + (rigth.getFy() - sol1y)
							* (rigth.getFy() - sol1y));
					double d2 = Math.sqrt((rigth.getFx() - sol2x)
							* (rigth.getFx() - sol2x) + (rigth.getFy() - sol2y)
							* (rigth.getFy() - sol2y));

					if (d1 < d2) {

						if (!left.getFlag()) {
							left.setFx(sol2x);
							left.setFy(sol2y);
							if ((left.getType() == 6)) {
								left.setFlag(true);
							}
						}
					} else {

						if (!left.getFlag()) {
							left.setFx(sol1x);
							left.setFy(sol1y);
							if ((left.getType() == 6)) {
								left.setFlag(true);
							}
						}
					}

					double solx1 = 0;
					double soly1 = 0;
					double solx2 = 0;
					double soly2 = 0;

					double xo1 = ref.getFx();
					double yo1 = ref.getFy();
					double xo2 = rigth.getFx();
					double yo2 = rigth.getFy();

					double r1 = 1.2;
					double r2 = 1.7;

					double a1 = -2.0 * xo1;
					double b1 = -2.0 * yo1;
					double a2 = -2.0 * xo2;
					double b2 = -2.0 * yo2;
					double a12 = a1 - a2;
					double b12 = b1 - b2;
					double aa = a12 * a12 + b12 * b12;

					double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
					double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
					double C12 = C1 - C2;
					double bb, cc;
					boolean bpga = Math.abs(b12) >= Math.abs(a12);

					if (bpga) {
						bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
						cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
					} else {
						bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
						cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
					}
					double delta = bb * bb - 4.0 * aa * cc;

					solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
					solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta)))
							/ aa;

					if (bpga) {
						soly1 = -(a12 * solx1 + C12) / b12;
						soly2 = -(a12 * solx2 + C12) / b12;
					} else {
						soly1 = solx1;
						solx1 = -(b12 * solx1 + C12) / a12;
						soly2 = solx2;
						solx2 = -(b12 * solx2 + C12) / a12;
					}

					if (soly1 > soly2) {
						if (!top.getFlag()) {
							top.setFx(solx1);
							top.setFy(soly1);
							if ((top.getType() == 6)) {
								top.setFlag(true);
							}
						}
						if (!bottom.getFlag()) {
							bottom.setFx(solx2);
							bottom.setFy(soly2);
							if ((bottom.getType() == 6)) {
								bottom.setFlag(true);
							}
						}
					} else {
						if (!top.getFlag()) {
							top.setFx(solx2);
							top.setFy(soly2);
							if ((top.getType() == 6)) {
								top.setFlag(true);
							}
						}
						if (!bottom.getFlag()) {
							bottom.setFx(solx1);
							bottom.setFy(soly1);
							if ((bottom.getType() == 6)) {
								bottom.setFlag(true);
							}
						}

					}

				} else {

					double sol1x = 0;
					double sol1y = 0;
					double sol2x = 0;
					double sol2y = 0;

					double x1 = ref.getFx();
					double y1 = ref.getFy();
					double x2 = left.getFx();
					double y2 = left.getFy();
					double a = y2 - y1;
					double b = x1 - x2;
					double c = y1 * x2 - x1 * y2;

					double xo = x1;
					double yo = y1;
					double r = 1.2;

					if (b != 0.0) {
						double u = a * (c + b * yo) - b * b * xo;
						double v = a * a + b * b;
						double w = c + b * yo;
						double deltap = u * u - v
								* (b * b * (xo * xo - r * r) + w * w);
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

					double d1 = Math.sqrt((left.getFx() - sol1x)
							* (left.getFx() - sol1x) + (left.getFy() - sol1y)
							* (left.getFy() - sol1y));
					double d2 = Math.sqrt((left.getFx() - sol2x)
							* (left.getFx() - sol2x) + (left.getFy() - sol2y)
							* (left.getFy() - sol2y));

					if (d1 < d2) {

						if (!rigth.getFlag()) {
							rigth.setFx(sol2x);
							rigth.setFy(sol2y);
							if ((rigth.getType() == 6)) {
								rigth.setFlag(true);
							}
						}
					} else {

						if (!rigth.getFlag()) {
							rigth.setFx(sol1x);
							rigth.setFy(sol1y);
							if ((rigth.getType() == 6)) {
								rigth.setFlag(true);
							}
						}
					}

					double solx1 = 0;
					double soly1 = 0;
					double solx2 = 0;
					double soly2 = 0;

					double xo1 = ref.getFx();
					double yo1 = ref.getFy();
					double xo2 = left.getFx();
					double yo2 = left.getFy();

					double r1 = 1.2;
					double r2 = 1.7;

					double a1 = -2.0 * xo1;
					double b1 = -2.0 * yo1;
					double a2 = -2.0 * xo2;
					double b2 = -2.0 * yo2;
					double a12 = a1 - a2;
					double b12 = b1 - b2;
					double aa = a12 * a12 + b12 * b12;

					double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
					double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
					double C12 = C1 - C2;
					double bb, cc;
					boolean bpga = Math.abs(b12) >= Math.abs(a12);

					if (bpga) {
						bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
						cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
					} else {
						bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
						cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
					}
					double delta = bb * bb - 4.0 * aa * cc;

					solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
					solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta)))
							/ aa;

					if (bpga) {
						soly1 = -(a12 * solx1 + C12) / b12;
						soly2 = -(a12 * solx2 + C12) / b12;
					} else {
						soly1 = solx1;
						solx1 = -(b12 * solx1 + C12) / a12;
						soly2 = solx2;
						solx2 = -(b12 * solx2 + C12) / a12;
					}

					if (soly1 > soly2) {
						if (!top.getFlag()) {
							top.setFx(solx1);
							top.setFy(soly1);
							if ((top.getType() == 6)) {
								top.setFlag(true);
							}
						}
						if (!bottom.getFlag()) {
							bottom.setFx(solx2);
							bottom.setFy(soly2);
							if ((bottom.getType() == 6)) {
								bottom.setFlag(true);
							}
						}
					} else {
						if (!top.getFlag()) {
							top.setFx(solx2);
							top.setFy(soly2);
							if ((top.getType() == 6)) {
								top.setFlag(true);
							}
						}
						if (!bottom.getFlag()) {
							bottom.setFx(solx1);
							bottom.setFy(soly1);
							if ((bottom.getType() == 6)) {
								bottom.setFlag(true);
							}
						}

					}

				}
			}
		}
	}

	public void twoCercles3(Atom ref, Atom[] tab) {

		Atom top = tab[0];
		Atom rigth = tab[1];
		Atom bottom = tab[2];
		Atom left = tab[3];

		if (top.getFlag()) {
			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = top.getFx();
			double y2 = top.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double d1 = Math.sqrt((top.getFx() - sol1x) * (top.getFx() - sol1x)
					+ (top.getFy() - sol1y) * (top.getFy() - sol1y));
			double d2 = Math.sqrt((top.getFx() - sol2x) * (top.getFx() - sol2x)
					+ (top.getFy() - sol2y) * (top.getFy() - sol2y));

			if (d1 < d2) {
				if (!top.getFlag()) {
					top.setFx(sol1x);
					top.setFy(sol1y);
					if ((top.getType() == 6)) {
						top.setFlag(true);
					}
				}
				if (!bottom.getFlag()) {
					bottom.setFx(sol2x);
					bottom.setFy(sol2y);
					if ((bottom.getType() == 6)) {
						bottom.setFlag(true);
					}
				}
			} else {
				if (!top.getFlag()) {
					top.setFx(sol2x);
					top.setFy(sol2y);
					if ((top.getType() == 6)) {
						top.setFlag(true);
					}
				}
				if (!bottom.getFlag()) {
					bottom.setFx(sol1x);
					bottom.setFy(sol1y);
					if ((bottom.getType() == 6)) {
						bottom.setFlag(true);
					}
				}
			}

			double solx1 = 0;
			double soly1 = 0;
			double solx2 = 0;
			double soly2 = 0;

			double xo1 = ref.getFx();
			double yo1 = ref.getFy();
			double xo2 = top.getFx();
			double yo2 = top.getFy();

			double r1 = 1.2;
			double r2 = 2;

			double a1 = -2.0 * xo1;
			double b1 = -2.0 * yo1;
			double a2 = -2.0 * xo2;
			double b2 = -2.0 * yo2;
			double a12 = a1 - a2;
			double b12 = b1 - b2;
			double aa = a12 * a12 + b12 * b12;

			double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
			double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
			double C12 = C1 - C2;
			double bb, cc;
			boolean bpga = Math.abs(b12) >= Math.abs(a12);

			if (bpga) {
				bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
				cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
			} else {
				bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
				cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
			}
			double delta = bb * bb - 4.0 * aa * cc;

			solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
			solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta))) / aa;

			if (bpga) {
				soly1 = -(a12 * solx1 + C12) / b12;
				soly2 = -(a12 * solx2 + C12) / b12;
			} else {
				soly1 = solx1;
				solx1 = -(b12 * solx1 + C12) / a12;
				soly2 = solx2;
				solx2 = -(b12 * solx2 + C12) / a12;
			}

			if (solx1 > solx2) {
				if (!rigth.getFlag()) {
					rigth.setFx(solx1);
					rigth.setFy(soly1);
					if ((rigth.getType() == 6)) {
						rigth.setFlag(true);
					}
				}
				if (!left.getFlag()) {
					left.setFx(solx2);
					left.setFy(soly2);
					if ((left.getType() == 6)) {
						left.setFlag(true);
					}
				}
			} else {
				if (!rigth.getFlag()) {
					rigth.setFx(solx2);
					rigth.setFy(soly2);
					if ((rigth.getType() == 6)) {
						rigth.setFlag(true);
					}
				}
				if (!left.getFlag()) {
					left.setFx(solx1);
					left.setFy(soly1);
					if ((left.getType() == 6)) {
						left.setFlag(true);
					}
				}

			}

		} else {
			if (bottom.getFlag()) {

				double sol1x = 0;
				double sol1y = 0;
				double sol2x = 0;
				double sol2y = 0;

				double x1 = ref.getFx();
				double y1 = ref.getFy();
				double x2 = bottom.getFx();
				double y2 = bottom.getFy();
				double a = y2 - y1;
				double b = x1 - x2;
				double c = y1 * x2 - x1 * y2;

				double xo = x1;
				double yo = y1;
				double r = 1.2;

				if (b != 0.0) {
					double u = a * (c + b * yo) - b * b * xo;
					double v = a * a + b * b;
					double w = c + b * yo;
					double deltap = u * u - v
							* (b * b * (xo * xo - r * r) + w * w);
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

				double d1 = Math.sqrt((bottom.getFx() - sol1x)
						* (bottom.getFx() - sol1x) + (bottom.getFy() - sol1y)
						* (bottom.getFy() - sol1y));
				double d2 = Math.sqrt((bottom.getFx() - sol2x)
						* (bottom.getFx() - sol2x) + (bottom.getFy() - sol2y)
						* (bottom.getFy() - sol2y));

				if (d1 < d2) {

					if (!top.getFlag()) {
						top.setFx(sol2x);
						top.setFy(sol2y);
						if ((top.getType() == 6)) {
							top.setFlag(true);
						}
					}
				} else {

					if (!top.getFlag()) {
						top.setFx(sol1x);
						top.setFy(sol1y);
						if ((top.getType() == 6)) {
							top.setFlag(true);
						}
					}
				}

				double solx1 = 0;
				double soly1 = 0;
				double solx2 = 0;
				double soly2 = 0;

				double xo1 = ref.getFx();
				double yo1 = ref.getFy();
				double xo2 = bottom.getFx();
				double yo2 = bottom.getFy();

				double r1 = 1.2;
				double r2 = 2;

				double a1 = -2.0 * xo1;
				double b1 = -2.0 * yo1;
				double a2 = -2.0 * xo2;
				double b2 = -2.0 * yo2;
				double a12 = a1 - a2;
				double b12 = b1 - b2;
				double aa = a12 * a12 + b12 * b12;

				double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
				double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
				double C12 = C1 - C2;
				double bb, cc;
				boolean bpga = Math.abs(b12) >= Math.abs(a12);

				if (bpga) {
					bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
					cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
				} else {
					bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
					cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
				}
				double delta = bb * bb - 4.0 * aa * cc;

				solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
				solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta))) / aa;

				if (bpga) {
					soly1 = -(a12 * solx1 + C12) / b12;
					soly2 = -(a12 * solx2 + C12) / b12;
				} else {
					soly1 = solx1;
					solx1 = -(b12 * solx1 + C12) / a12;
					soly2 = solx2;
					solx2 = -(b12 * solx2 + C12) / a12;
				}

				if (solx1 > solx2) {
					if (!rigth.getFlag()) {
						rigth.setFx(solx1);
						rigth.setFy(soly1);
						if ((rigth.getType() == 6)) {
							rigth.setFlag(true);
						}
					}
					if (!left.getFlag()) {
						left.setFx(solx2);
						left.setFy(soly2);
						if ((left.getType() == 6)) {
							left.setFlag(true);
						}
					}
				} else {
					if (!rigth.getFlag()) {
						rigth.setFx(solx2);
						rigth.setFy(soly2);
						if ((rigth.getType() == 6)) {
							rigth.setFlag(true);
						}
					}
					if (!left.getFlag()) {
						left.setFx(solx1);
						left.setFy(soly1);
						if ((left.getType() == 6)) {
							left.setFlag(true);
						}
					}

				}

			} else {
				if (rigth.getFlag()) {

					double sol1x = 0;
					double sol1y = 0;
					double sol2x = 0;
					double sol2y = 0;

					double x1 = ref.getFx();
					double y1 = ref.getFy();
					double x2 = rigth.getFx();
					double y2 = rigth.getFy();
					double a = y2 - y1;
					double b = x1 - x2;
					double c = y1 * x2 - x1 * y2;

					double xo = x1;
					double yo = y1;
					double r = 1.2;

					if (b != 0.0) {
						double u = a * (c + b * yo) - b * b * xo;
						double v = a * a + b * b;
						double w = c + b * yo;
						double deltap = u * u - v
								* (b * b * (xo * xo - r * r) + w * w);
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

					double d1 = Math.sqrt((rigth.getFx() - sol1x)
							* (rigth.getFx() - sol1x) + (rigth.getFy() - sol1y)
							* (rigth.getFy() - sol1y));
					double d2 = Math.sqrt((rigth.getFx() - sol2x)
							* (rigth.getFx() - sol2x) + (rigth.getFy() - sol2y)
							* (rigth.getFy() - sol2y));

					if (d1 < d2) {

						if (!left.getFlag()) {
							left.setFx(sol2x);
							left.setFy(sol2y);
							if ((left.getType() == 6)) {
								left.setFlag(true);
							}
						}
					} else {

						if (!left.getFlag()) {
							left.setFx(sol1x);
							left.setFy(sol1y);
							if ((left.getType() == 6)) {
								left.setFlag(true);
							}
						}
					}

					double solx1 = 0;
					double soly1 = 0;
					double solx2 = 0;
					double soly2 = 0;

					double xo1 = ref.getFx();
					double yo1 = ref.getFy();
					double xo2 = rigth.getFx();
					double yo2 = rigth.getFy();

					double r1 = 1.2;
					double r2 = 2;

					double a1 = -2.0 * xo1;
					double b1 = -2.0 * yo1;
					double a2 = -2.0 * xo2;
					double b2 = -2.0 * yo2;
					double a12 = a1 - a2;
					double b12 = b1 - b2;
					double aa = a12 * a12 + b12 * b12;

					double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
					double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
					double C12 = C1 - C2;
					double bb, cc;
					boolean bpga = Math.abs(b12) >= Math.abs(a12);

					if (bpga) {
						bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
						cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
					} else {
						bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
						cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
					}
					double delta = bb * bb - 4.0 * aa * cc;

					solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
					solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta)))
							/ aa;

					if (bpga) {
						soly1 = -(a12 * solx1 + C12) / b12;
						soly2 = -(a12 * solx2 + C12) / b12;
					} else {
						soly1 = solx1;
						solx1 = -(b12 * solx1 + C12) / a12;
						soly2 = solx2;
						solx2 = -(b12 * solx2 + C12) / a12;
					}

					if (soly1 > soly2) {
						if (!top.getFlag()) {
							top.setFx(solx1);
							top.setFy(soly1);
							if ((top.getType() == 6)) {
								top.setFlag(true);
							}
						}
						if (!bottom.getFlag()) {
							bottom.setFx(solx2);
							bottom.setFy(soly2);
							if ((bottom.getType() == 6)) {
								bottom.setFlag(true);
							}
						}
					} else {
						if (!top.getFlag()) {
							top.setFx(solx2);
							top.setFy(soly2);
							if ((top.getType() == 6)) {
								top.setFlag(true);
							}
						}
						if (!bottom.getFlag()) {
							bottom.setFx(solx1);
							bottom.setFy(soly1);
							if ((bottom.getType() == 6)) {
								bottom.setFlag(true);
							}
						}

					}

				} else {

					double sol1x = 0;
					double sol1y = 0;
					double sol2x = 0;
					double sol2y = 0;

					double x1 = ref.getFx();
					double y1 = ref.getFy();
					double x2 = left.getFx();
					double y2 = left.getFy();
					double a = y2 - y1;
					double b = x1 - x2;
					double c = y1 * x2 - x1 * y2;

					double xo = x1;
					double yo = y1;
					double r = 1.2;

					if (b != 0.0) {
						double u = a * (c + b * yo) - b * b * xo;
						double v = a * a + b * b;
						double w = c + b * yo;
						double deltap = u * u - v
								* (b * b * (xo * xo - r * r) + w * w);
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

					double d1 = Math.sqrt((left.getFx() - sol1x)
							* (left.getFx() - sol1x) + (left.getFy() - sol1y)
							* (left.getFy() - sol1y));
					double d2 = Math.sqrt((left.getFx() - sol2x)
							* (left.getFx() - sol2x) + (left.getFy() - sol2y)
							* (left.getFy() - sol2y));

					if (d1 < d2) {

						if (!rigth.getFlag()) {
							rigth.setFx(sol2x);
							rigth.setFy(sol2y);
							if ((rigth.getType() == 6)) {
								rigth.setFlag(true);
							}
						}
					} else {

						if (!rigth.getFlag()) {
							rigth.setFx(sol1x);
							rigth.setFy(sol1y);
							if ((rigth.getType() == 6)) {
								rigth.setFlag(true);
							}
						}
					}

					double solx1 = 0;
					double soly1 = 0;
					double solx2 = 0;
					double soly2 = 0;

					double xo1 = ref.getFx();
					double yo1 = ref.getFy();
					double xo2 = left.getFx();
					double yo2 = left.getFy();

					double r1 = 1.2;
					double r2 = 2;

					double a1 = -2.0 * xo1;
					double b1 = -2.0 * yo1;
					double a2 = -2.0 * xo2;
					double b2 = -2.0 * yo2;
					double a12 = a1 - a2;
					double b12 = b1 - b2;
					double aa = a12 * a12 + b12 * b12;

					double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
					double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
					double C12 = C1 - C2;
					double bb, cc;
					boolean bpga = Math.abs(b12) >= Math.abs(a12);

					if (bpga) {
						bb = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
						cc = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
					} else {
						bb = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
						cc = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
					}
					double delta = bb * bb - 4.0 * aa * cc;

					solx1 = -0.5 * (bb + 1.0 * Math.sqrt(Math.abs(delta))) / aa;
					solx2 = -0.5 * (bb + -1.0 * Math.sqrt(Math.abs(delta)))
							/ aa;

					if (bpga) {
						soly1 = -(a12 * solx1 + C12) / b12;
						soly2 = -(a12 * solx2 + C12) / b12;
					} else {
						soly1 = solx1;
						solx1 = -(b12 * solx1 + C12) / a12;
						soly2 = solx2;
						solx2 = -(b12 * solx2 + C12) / a12;
					}

					if (soly1 > soly2) {
						if (!top.getFlag()) {
							top.setFx(solx1);
							top.setFy(soly1);
							if ((top.getType() == 6)) {
								top.setFlag(true);
							}
						}
						if (!bottom.getFlag()) {
							bottom.setFx(solx2);
							bottom.setFy(soly2);
							if ((bottom.getType() == 6)) {
								bottom.setFlag(true);
							}
						}
					} else {
						if (!top.getFlag()) {
							top.setFx(solx2);
							top.setFy(soly2);
							if ((top.getType() == 6)) {
								top.setFlag(true);
							}
						}
						if (!bottom.getFlag()) {
							bottom.setFx(solx1);
							bottom.setFy(soly1);
							if ((bottom.getType() == 6)) {
								bottom.setFlag(true);
							}
						}

					}

				}
			}
		}
	}

	public LinkedList interDC(double xo, double yo, double r, double x, double y) {

		LinkedList res = new LinkedList();

		double a = y - yo;
		double b = xo - x;
		double c = yo * x - xo * y;
		double sol1x = 0;
		double sol2x = 0;
		double sol1y = 0;
		double sol2y = 0;

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

		double d1 = Math.sqrt((x - sol1x) * (x - sol1x) + (y - sol1y)
				* (y - sol1y));
		double d2 = Math.sqrt((x - sol2x) * (x - sol2x) + (y - sol2y)
				* (y - sol2y));

		if (d1 < d2) {
			res.add(sol1x);
			res.add(sol1y);
		} else {
			res.add(sol2x);
			res.add(sol2y);
		}

		return res;

	}

	public void twoCercles(Atom fixed, Atom ref, Atom at1, Atom at2, Atom at3) {

		if (!fixed.getFlag()) {
			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x1 = ref.getFx();
			double y1 = ref.getFy();
			double x2 = fixed.getFx();
			double y2 = fixed.getFy();
			double a = y2 - y1;
			double b = x1 - x2;
			double c = y1 * x2 - x1 * y2;

			double xo = x1;
			double yo = y1;
			double r = 1.2;

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

			double d1 = Math.sqrt((fixed.getFx() - sol1x)
					* (fixed.getFx() - sol1x) + (fixed.getFy() - sol1y)
					* (fixed.getFy() - sol1y));
			double d2 = Math.sqrt((fixed.getFx() - sol2x)
					* (fixed.getFx() - sol2x) + (fixed.getFy() - sol2y)
					* (fixed.getFy() - sol2y));

			if (d1 < d2) {
				fixed.setFx(sol1x);
				fixed.setFy(sol1y);
				if ((fixed.getType() == 5) || (fixed.getType() == 6)) {
					fixed.setFlag(true);
				}
			} else {
				fixed.setFx(sol2x);
				fixed.setFy(sol2y);
				if ((fixed.getType() == 5) || (fixed.getType() == 6)) {
					fixed.setFlag(true);
				}
			}

		}

		double xo1 = ref.getFx();
		double yo1 = ref.getFy();
		double xo2 = fixed.getFx();
		double yo2 = fixed.getFy();

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		double r1 = 1.2;
		double r2 = 1.7;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo2;
		double b2 = -2.0 * yo2;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double a = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
		double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
		double C12 = C1 - C2;
		double b, c;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);

		if (bpga) {
			b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = b * b - 4.0 * a * c;

		solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
		solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			solx1 = -(b12 * solx1 + C12) / a12;
			soly2 = solx2;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		double solx;
		double soly;

		double d1 = Math.sqrt((at1.getFx() - solx1) * (at1.getFx() - solx1)
				+ (at1.getFy() - soly1) * (at1.getFy() - soly1));
		double d2 = Math.sqrt((at2.getFx() - solx1) * (at2.getFx() - solx1)
				+ (at2.getFy() - soly1) * (at2.getFy() - soly1));
		double d3 = Math.sqrt((at3.getFx() - solx1) * (at3.getFx() - solx1)
				+ (at3.getFy() - soly1) * (at3.getFy() - soly1));

		LinkedList trait = new LinkedList();

		if ((d1 < d2) && (d1 < d3)) {
			if (!at1.getFlag()) {
				at1.setFx(solx1);
				at1.setFy(soly1);
				trait.add(at2);
				trait.add(at3);
				if ((at1.getType() == 5) || (at1.getType() == 6)) {
					at1.setFlag(true);
				}
			}
		} else {
			if ((d2 < d1) && (d2 < d3)) {
				if (!at2.getFlag()) {
					at2.setFx(solx1);
					at2.setFy(soly1);
					trait.add(at1);
					trait.add(at3);
					if ((at3.getType() == 5) || (at3.getType() == 6)) {
						at3.setFlag(true);
					}
				}
			} else {
				if (!at3.getFlag()) {
					at3.setFx(solx1);
					at3.setFy(soly1);
					trait.add(at1);
					trait.add(at2);
					if ((at3.getType() == 5) || (at3.getType() == 6)) {
						at3.setFlag(true);
					}
				}
			}
		}

		Atom att1 = (Atom) trait.get(0);

		Atom att2 = (Atom) trait.get(0);
		Atom last = null;

		double d4 = Math.sqrt((att1.getFx() - solx2) * (att1.getFx() - solx2)
				+ (att1.getFy() - soly2) * (att1.getFy() - soly2));
		double d5 = Math.sqrt((att2.getFx() - solx2) * (att2.getFx() - solx2)
				+ (att2.getFy() - soly2) * (att2.getFy() - soly2));

		if (d4 < d5) {
			if (!att1.getFlag()) {
				att1.setFx(solx2);
				att1.setFy(soly2);
				last = att2;
				if ((att1.getType() == 5) || (att1.getType() == 6)) {
					att1.setFlag(true);
				}
			}

		} else {
			if (!att2.getFlag()) {
				att2.setFx(solx2);
				att2.setFy(soly2);
				last = att1;
				if ((att2.getType() == 5) || (att2.getType() == 6)) {
					att2.setFlag(true);
				}
			}

		}

		Atom s1;
		Atom s2;
		if (last == at1) {
			s1 = at2;
			s2 = at3;
		} else {
			if (last == at2) {
				s1 = at1;
				s2 = at3;
			} else {
				s1 = at1;
				s2 = at3;
			}
		}

		double aa = s1.getFy() - s2.getFy();
		double bb = s2.getFx() - s1.getFx();
		double cc = -(a * s1.getFx() + b * s1.getFy());
		double u = (aa * last.getFx() + bb * last.getFy() + cc)
				/ (aa * aa + bb * bb);
		double x = last.getFx() - u * aa;
		double y = last.getFy() - u * bb;

		x = 2.0 * x - fixed.getFx();
		y = 2.0 * y - fixed.getFy();

		if (!last.getFlag()) {
			last.setFx(x);
			last.setFy(y);
			if ((last.getType() == 5) || (last.getType() == 6)) {
				last.setFlag(true);
			}
		}

	}

	public void twoCercles2(Atom ref, Atom fixed, Atom at1, Atom at2) {

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		double xo1 = ref.getFx();
		double yo1 = ref.getFy();
		double xo2 = fixed.getFx();
		double yo2 = fixed.getFy();

		double r1 = 1.2;
		double r2 = 2.43;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo2;
		double b2 = -2.0 * yo2;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double a = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
		double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
		double C12 = C1 - C2;
		double b, c;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);

		if (bpga) {
			b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = b * b - 4.0 * a * c;

		solx1 = -0.5 * (b + 1.0 * Math.sqrt(Math.abs(delta))) / a;
		solx2 = -0.5 * (b + -1.0 * Math.sqrt(Math.abs(delta))) / a;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			solx1 = -(b12 * solx1 + C12) / a12;
			soly2 = solx2;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		double angle1 = angleOld(at1, ref, at2);
		double angle2 = angleNew(solx1, soly1, ref, solx2, soly2);

		if (((angle1 < 0) && (angle2 < 0)) || ((angle1 > 0) && (angle2 > 0))) {
			if (!at1.getFlag()) {
				at1.setFx(solx1);
				at1.setFy(soly1);
				if ((at1.getType() == 5) || (at1.getType() == 6)) {
					at1.setFlag(true);
				}
			}
			if (!at2.getFlag()) {
				at2.setFx(solx2);
				at2.setFy(soly2);
				if ((at2.getType() == 5) || (at2.getType() == 6)) {
					at2.setFlag(true);
				}
			}
		} else {
			if (!at1.getFlag()) {
				at1.setFx(solx2);
				at1.setFy(soly2);
				if ((at1.getType() == 5) || (at1.getType() == 6)) {
					at1.setFlag(true);
				}
			}
			if (!at2.getFlag()) {
				at2.setFx(solx1);
				at2.setFy(soly1);
				if ((at2.getType() == 5) || (at2.getType() == 6)) {
					at2.setFlag(true);
				}
			}
		}

	}

	public void twoCercles2(Atom ref, Atom fixed, Atom at1, Atom at2, Atom at3) {

		Atom[] tab = new Atom[4];

		tab[0] = fixed;
		tab[1] = at1;
		tab[2] = at2;
		tab[3] = at3;

		Atom[] tab2 = rigthOrder(ref, tab);

		twoCercles3(ref, tab2);

	}

	public Atom getPredecesor(Atom at) {

		LinkedList neig = at.getNeig();
		if (neig.size() == 0)
			return null;
		for (int i = 0; i < neig.size(); i++) {
			Atom a = (Atom) neig.get(i);
			if ((!a.getTerm()) && (a.getFlag()))
				return a;
		}
		return null;
	}

	public Atom getPredecesor(Ring ring) {

		LinkedList at = ring.getAtoms();

		for (int i = 0; i < at.size(); i++) {
			Atom a = getAtom((Integer) at.get(i));
			LinkedList neig = a.getNeig();
			for (int j = 0; j < neig.size(); j++) {
				Atom aa = (Atom) neig.get(j);
				if ((!aa.getTerm()) && (aa.getFlag()))
					return aa;
			}
		}
		return null;
	}

	public void traitChaine(Atom deb, Chaine chaine) {

		LinkedList atomId = chaine.getAtoms();
		int order = 0;
		if ((Integer) atomId.get(0) == deb.getId()) {
			order = 1;
		} else {
			if ((Integer) atomId.get(atomId.size() - 1) == deb.getId()) {
				order = -1;
			}
		}

		if (order == 0)
			return;
		if (order == 1) {

			Atom first = getAtom((Integer) atomId.get(0));
			Atom second = getAtom((Integer) atomId.get(1));
			Atom last;

			Atom prede = getPredecesor(first);
			Atom prede2 = getPredecesor(getAtom((Integer) atomId.get(atomId
					.size() - 1)));
			if (prede2 != null) {

				traitChaine(getAtom((Integer) atomId.get(atomId.size() - 1)),
						chaine);
			} else {

				if (prede != null) {
					last = prede.clone();
					double nx = 2.0 * first.getFx() - prede.getFx();
					double ny = 2.0 * first.getFy() - prede.getFy();
					last.setFx(nx);
					last.setFy(ny);

				} else {

					if (atomId.size() % 2 == 0) {
						last = getAtom((Integer) atomId.get(atomId.size() - 2));
					} else {
						last = getAtom((Integer) atomId.get(atomId.size() - 1));
					}
				}

				int sens;

				Atom big = getBig(chaine);

				int ind = 0;

				for (int i = 1; i < atomId.size() - 1; i++) {
					int idTemp = (Integer) atomId.get(i);
					if (idTemp == big.getId()) {
						ind = i;
						continue;
					}

				}

				int type = 0;

				double angle;

				if ((big.getId() == first.getId())
						|| (big.getId() == last.getId())) {
					angle = angle = angle(last, first, second);
				} else {
					angle = angle(last, first, big);
				}

				if (angle >= 0) {
					sens = 1;
				} else {
					sens = 0;
				}

				if (ind % 2 != 0) {
					if (sens == 1)
						sens = 0;
					else
						sens = 1;
				}

				// pairs
				int cont = 1;
				for (int i = 2; i < atomId.size(); i = i + 2) {

					Atom trait = getAtom((Integer) atomId.get(i));

					intersection(first, last, cont, trait);
					trait.setFlag(true);

					if (i != atomId.size() - 1) {
						LinkedList neig = trait.getNeig();
						LinkedList temp = new LinkedList();
						LinkedList temp2 = new LinkedList();
						for (int s = 0; s < neig.size(); s++) {
							Atom aa = (Atom) neig.get(s);

							if (atomId.indexOf(aa.getId()) == -1) {
								temp.add(aa);
							} else
								temp2.add(aa);
						}
						if (temp.size() == 1) {
							Atom att = (Atom) temp.get(0);
							if (!att.getFlag()) {
								double x0 = trait.getFx();
								double y0 = trait.getFy();
								double x1 = first.getFx();
								double y1 = first.getFy();
								double x2 = last.getFx();
								double y2 = last.getFy();
								double a = y2 - y1;
								double b = x1 - x2;
								double c = y1 * x2 - x1 * y2;
								double a2 = -b;
								double b2 = a;
								double c2 = b * x0 - a * y0;
								traitAtomChaine(att, x0, y0, a2, b2, c2, sens);
								traitAtom(att);
								if ((att.getType() == 5)
										|| (att.getType() == 6)) {
									att.setFlag(true);
								}
							}

						}

						if (temp.size() == 2) {

							Atom att = (Atom) temp.get(0);
							Atom tempo = att.clone();
							Atom att2 = (Atom) temp.get(1);

							double x0 = trait.getFx();
							double y0 = trait.getFy();
							double x1 = first.getFx();
							double y1 = first.getFy();
							double x2 = last.getFx();
							double y2 = last.getFy();
							double a = y2 - y1;
							double b = x1 - x2;
							double c = y1 * x2 - x1 * y2;
							double a2 = -b;
							double b2 = a;
							double c2 = b * x0 - a * y0;
							traitAtomChaine(tempo, x0, y0, a2, b2, c2, sens);

							double OAx = tempo.getFx() - trait.getFx();
							double OAy = tempo.getFy() - trait.getFy();

							double resx = trait.getFx() + Math.cos(0.69) * OAx
									- Math.sin(0.69) * OAy;
							double resy = trait.getFy() + Math.sin(0.69) * OAx
									+ Math.cos(0.69) * OAy;

							double resx2 = trait.getFx() + Math.cos(-0.69)
									* OAx - Math.sin(-0.69) * OAy;
							double resy2 = trait.getFy() + Math.sin(-0.69)
									* OAx + Math.cos(-0.69) * OAy;

							Atom ato1 = (Atom) temp2.get(0);
							Atom ato2 = (Atom) temp2.get(1);

							double angle1 = angleOld2(att, trait, att2);
							double angle2 = angleNew(resx, resy, trait, resx2,
									resy2);

							if ((angle1 < 0) && (angle2 < 0) || (angle1 > 0)
									&& (angle2 > 0)) {
								if (!att.getFlag()) {
									att.setFx(resx);
									att.setFy(resy);
									if ((att.getType() == 5)
											|| (att.getType() == 6)) {
										att.setFlag(true);
									}
								}
								if (!att2.getFlag()) {
									att2.setFx(resx2);
									att2.setFy(resy2);
									if ((att2.getType() == 5)
											|| (att2.getType() == 6)) {
										att2.setFlag(true);
									}
								}

							} else {
								if (!att.getFlag()) {
									att.setFx(resx2);
									att.setFy(resy2);
									if ((att.getType() == 5)
											|| (att.getType() == 6)) {
										att.setFlag(true);
									}
								}
								if (!att2.getFlag()) {
									att2.setFx(resx);
									att2.setFy(resy);
									if ((att2.getType() == 5)
											|| (att2.getType() == 6)) {
										att2.setFlag(true);
									}
								}

							}

						}

					}

					cont = cont + 1;

				}

				// impairs
				for (int i = 1; i < atomId.size(); i = i + 2) {
					Atom trait = getAtom((Integer) atomId.get(i));
					if (i == 1) {

						Atom pred = getAtom((Integer) atomId.get(i - 1));
						Atom suc = getAtom((Integer) atomId.get(i + 1));
						intersectionImpair(trait, pred, suc, sens);
						trait.setFlag(true);

						double angletype = angle(first, trait, last);
						if (angletype > 0)
							type = -1;
						else
							type = 1;
					} else {

						Atom pred = getAtom((Integer) atomId.get(i - 1));
						Atom pred2 = getAtom((Integer) atomId.get(i - 2));
						intersectionImpair2(trait, pred, pred2, type);
						trait.setFlag(true);

					}

					if (i != atomId.size() - 1) {
						LinkedList neig = trait.getNeig();
						LinkedList temp = new LinkedList();
						LinkedList temp2 = new LinkedList();
						for (int s = 0; s < neig.size(); s++) {
							Atom aa = (Atom) neig.get(s);
							if (atomId.indexOf(aa.getId()) == -1) {
								temp.add(aa);
							} else
								temp2.add(aa);
						}
						if (temp.size() == 1) {
							Atom att = (Atom) temp.get(0);
							if (!att.getFlag()) {
								double x0 = trait.getFx();
								double y0 = trait.getFy();
								double x1 = first.getFx();
								double y1 = first.getFy();
								double x2 = last.getFx();
								double y2 = last.getFy();
								double a = y2 - y1;
								double b = x1 - x2;
								double c = y1 * x2 - x1 * y2;
								double a2 = -b;
								double b2 = a;
								double c2 = b * x0 - a * y0;
								traitAtomChaine(att, x0, y0, a2, b2, c2, a, b,
										c);
								traitAtom(att);
								if ((att.getType() == 5)
										|| (att.getType() == 6)) {
									att.setFlag(true);
								}

							}
						}
						if (temp.size() == 2) {

							Atom att = (Atom) temp.get(0);
							Atom tempo = att.clone();
							Atom att2 = (Atom) temp.get(1);

							double x0 = trait.getFx();
							double y0 = trait.getFy();
							double x1 = first.getFx();
							double y1 = first.getFy();
							double x2 = last.getFx();
							double y2 = last.getFy();
							double a = y2 - y1;
							double b = x1 - x2;
							double c = y1 * x2 - x1 * y2;
							double a2 = -b;
							double b2 = a;
							double c2 = b * x0 - a * y0;
							traitAtomChaine(tempo, x0, y0, a2, b2, c2, a, b, c);

							double OAx = tempo.getFx() - trait.getFx();
							double OAy = tempo.getFy() - trait.getFy();

							double resx = trait.getFx() + Math.cos(0.69) * OAx
									- Math.sin(0.69) * OAy;
							double resy = trait.getFy() + Math.sin(0.69) * OAx
									+ Math.cos(0.69) * OAy;

							double resx2 = trait.getFx() + Math.cos(-0.69)
									* OAx - Math.sin(-0.69) * OAy;
							double resy2 = trait.getFy() + Math.sin(-0.69)
									* OAx + Math.cos(-0.69) * OAy;

							Atom ato1 = (Atom) temp2.get(0);
							Atom ato2 = (Atom) temp2.get(1);

							double angle1 = angleOld2(att, trait, att2);
							double angle2 = angleNew(resx, resy, trait, resx2,
									resy2);

							if ((angle1 < 0) && (angle2 < 0) || (angle1 > 0)
									&& (angle2 > 0)) {
								if (!att.getFlag()) {
									att.setFx(resx);
									att.setFy(resy);
									if ((att.getType() == 5)
											|| (att.getType() == 6)) {
										att.setFlag(true);
									}
								}
								if (!att2.getFlag()) {
									att2.setFx(resx2);
									att2.setFy(resy2);
									if ((att2.getType() == 5)
											|| (att2.getType() == 6)) {
										att2.setFlag(true);
									}
								}

							} else {
								if (!att.getFlag()) {
									att.setFx(resx2);
									att.setFy(resy2);
									if ((att.getType() == 5)
											|| (att.getType() == 6)) {
										att.setFlag(true);
									}
								}
								if (!att2.getFlag()) {
									att2.setFx(resx);
									att2.setFy(resy);
									if ((att2.getType() == 5)
											|| (att2.getType() == 6)) {
										att2.setFlag(true);
									}
								}

							}

						}

					}

				}

				LinkedList nt = first.getNeig();
				LinkedList nt2 = new LinkedList();

				if (nt.size() != 0) {
					for (int i = 0; i < nt.size(); i++) {
						Atom aaa = (Atom) nt.get(i);
						if (atomId.indexOf(aaa.getId()) == -1) {
							if (aaa.getType() == 6) {
								aaa.setFlag(false);
							}
							nt2.add(aaa);
						}
					}
					if (nt2.size() != 0) {
						if ((nt2.size() == 2) || (nt2.size() == 3)) {
							traitAtomC(first, getAtom((Integer) atomId.get(1)));
						}
						if (nt2.size() == 1) {
							Atom trait = (Atom) nt2.get(0);
							intersectionImpair3(trait, first,
									getAtom((Integer) atomId.get(1)),
									getAtom((Integer) atomId.get(2)));
						}
					}
				}

				first.setFlag(true);

				Atom ult = getAtom((Integer) atomId.get(atomId.size() - 1));

				LinkedList nnt = ult.getNeig();
				LinkedList nnt2 = new LinkedList();

				if (nnt.size() != 0) {
					for (int i = 0; i < nnt.size(); i++) {
						Atom aaa = (Atom) nnt.get(i);
						if (atomId.indexOf(aaa.getId()) == -1) {
							if (aaa.getType() == 6) {
								aaa.setFlag(false);
							}
							nnt2.add(aaa);
						}
					}
					if (nnt2.size() != 0) {
						if ((nnt2.size() == 2) || (nnt2.size() == 3))
							traitAtomC(ult, getAtom((Integer) atomId.get(atomId
									.size() - 2)));
						if (nnt2.size() == 1) {
							Atom trait = (Atom) nnt2.get(0);
							intersectionImpair3(trait, ult,
									getAtom((Integer) atomId
											.get(atomId.size() - 2)),
									getAtom((Integer) atomId
											.get(atomId.size() - 3)));
						}
					}
				}

			}
		}

		if (order == -1) {
			Atom first = getAtom((Integer) atomId.get(atomId.size() - 1));

			Atom second = getAtom((Integer) atomId.get(atomId.size() - 2));
			Atom last;

			Atom prede = getPredecesor(first);

			Atom prede2 = getPredecesor(getAtom((Integer) atomId.get(0)));

			if (prede2 != null) {

				traitChaine(getAtom((Integer) atomId.get(0)), chaine);
			} else {

				if (prede != null) {

					last = prede.clone();

					double nx = 2.0 * first.getFx() - prede.getFx();
					double ny = 2.0 * first.getFy() - prede.getFy();
					last.setFx(nx);
					last.setFy(ny);

				} else {

					if (atomId.size() % 2 == 0) {
						last = getAtom((Integer) atomId.get(1));
					} else {
						last = getAtom((Integer) atomId.get(0));
					}
				}

				int sens;

				Atom big = getBig(chaine);

				int ind = 0;

				for (int i = 1; i < atomId.size() - 1; i++) {
					int idTemp = (Integer) atomId.get(i);
					if (idTemp == big.getId()) {
						ind = i;
						continue;
					}

				}

				int type = 0;

				double angle;
				if ((big.getId() == first.getId())
						|| (big.getId() == last.getId())) {
					angle = angle = angle(last, first, second);
				} else {
					angle = angle(last, first, big);
				}

				if (angle >= 0) {
					sens = 1;
				} else {
					sens = 0;
				}

				if (ind % 2 != 0) {
					if (sens == 1)
						sens = 0;
					else
						sens = 1;
				}

				// pairs
				int cont = 1;
				for (int i = atomId.size() - 3; i >= 0; i = i - 2) {

					Atom trait = getAtom((Integer) atomId.get(i));

					intersection(first, last, cont, trait);
					trait.setFlag(true);

					if (i != 0) {
						LinkedList neig = trait.getNeig();
						LinkedList temp = new LinkedList();
						LinkedList temp2 = new LinkedList();
						for (int s = 0; s < neig.size(); s++) {
							Atom aa = (Atom) neig.get(s);
							if (atomId.indexOf(aa.getId()) == -1) {
								temp.add(aa);
							} else {
								temp2.add(aa);
							}
						}

						if (temp.size() == 1) {
							Atom att = (Atom) temp.get(0);
							if (!att.getFlag()) {
								double x0 = trait.getFx();
								double y0 = trait.getFy();
								double x1 = first.getFx();
								double y1 = first.getFy();
								double x2 = last.getFx();
								double y2 = last.getFy();
								double a = y2 - y1;
								double b = x1 - x2;
								double c = y1 * x2 - x1 * y2;
								double a2 = -b;
								double b2 = a;
								double c2 = b * x0 - a * y0;
								traitAtomChaine(att, x0, y0, a2, b2, c2, sens);
								traitAtom(att);
								if ((att.getType() == 5)
										|| (att.getType() == 6)) {
									att.setFlag(true);
								}

							}
						}

						if ((temp.size() == 2) && (temp2.size() == 2)) {

							Atom att = (Atom) temp.get(0);
							Atom tempo = att.clone();
							Atom att2 = (Atom) temp.get(1);

							double x0 = trait.getFx();
							double y0 = trait.getFy();
							double x1 = first.getFx();
							double y1 = first.getFy();
							double x2 = last.getFx();
							double y2 = last.getFy();
							double a = y2 - y1;
							double b = x1 - x2;
							double c = y1 * x2 - x1 * y2;
							double a2 = -b;
							double b2 = a;
							double c2 = b * x0 - a * y0;
							traitAtomChaine(tempo, x0, y0, a2, b2, c2, sens);

							double OAx = tempo.getFx() - trait.getFx();
							double OAy = tempo.getFy() - trait.getFy();

							double resx = trait.getFx() + Math.cos(0.69) * OAx
									- Math.sin(0.69) * OAy;
							double resy = trait.getFy() + Math.sin(0.69) * OAx
									+ Math.cos(0.69) * OAy;

							double resx2 = trait.getFx() + Math.cos(-0.69)
									* OAx - Math.sin(-0.69) * OAy;
							double resy2 = trait.getFy() + Math.sin(-0.69)
									* OAx + Math.cos(-0.69) * OAy;

							Atom ato1 = (Atom) temp2.get(0);
							Atom ato2 = (Atom) temp2.get(1);

							double angle1 = angleOld2(att, trait, att2);
							double angle2 = angleNew(resx, resy, trait, resx2,
									resy2);

							if ((angle1 < 0) && (angle2 < 0) || (angle1 > 0)
									&& (angle2 > 0)) {
								if (!att.getFlag()) {
									att.setFx(resx);
									att.setFy(resy);
									if ((att.getType() == 5)
											|| (att.getType() == 6)) {
										att.setFlag(true);
									}
								}
								if (!att2.getFlag()) {
									att2.setFx(resx2);
									att2.setFy(resy2);
									if ((att2.getType() == 5)
											|| (att2.getType() == 6)) {
										att2.setFlag(true);
									}
								}

							} else {
								if (!att.getFlag()) {
									att.setFx(resx2);
									att.setFy(resy2);
									if ((att.getType() == 5)
											|| (att.getType() == 6)) {
										att.setFlag(true);
									}
								}
								if (!att2.getFlag()) {
									att2.setFx(resx);
									att2.setFy(resy);
									if ((att2.getType() == 5)
											|| (att2.getType() == 6)) {
										att2.setFlag(true);
									}
								}

							}

						}

					}

					cont = cont + 1;

				}

				// impairs
				for (int i = atomId.size() - 2; i >= 0; i = i - 2) {

					Atom trait = getAtom((Integer) atomId.get(i));
					if (i == atomId.size() - 2) {
						Atom pred = getAtom((Integer) atomId.get(i - 1));
						Atom suc = getAtom((Integer) atomId.get(i + 1));
						intersectionImpair(trait, pred, suc, sens);
						trait.setFlag(true);

						double angletype = angle(first, trait, last);
						if (angletype > 0)
							type = -1;
						else
							type = 1;

					} else {
						Atom pred = getAtom((Integer) atomId.get(i + 1));
						Atom pred2 = getAtom((Integer) atomId.get(i + 2));
						intersectionImpair2(trait, pred, pred2, type);
						trait.setFlag(true);
					}

					if (i != 0) {
						LinkedList neig = trait.getNeig();
						LinkedList temp = new LinkedList();
						LinkedList temp2 = new LinkedList();
						for (int s = 0; s < neig.size(); s++) {
							Atom aa = (Atom) neig.get(s);
							if (atomId.indexOf(aa.getId()) == -1) {
								temp.add(aa);
							} else
								temp2.add(aa);
						}
						if (temp.size() == 1) {
							Atom att = (Atom) temp.get(0);
							if (!att.getFlag()) {
								double x0 = trait.getFx();
								double y0 = trait.getFy();
								double x1 = first.getFx();
								double y1 = first.getFy();
								double x2 = last.getFx();
								double y2 = last.getFy();
								double a = y2 - y1;
								double b = x1 - x2;
								double c = y1 * x2 - x1 * y2;
								double a2 = -b;
								double b2 = a;
								double c2 = b * x0 - a * y0;

								traitAtomChaine(att, x0, y0, a2, b2, c2, a, b,
										c);
								traitAtom(att);
								if ((att.getType() == 5)
										|| (att.getType() == 6)) {
									att.setFlag(true);
								}

							}
						}
						if (temp.size() == 2) {

							Atom att = (Atom) temp.get(0);
							Atom tempo = att.clone();
							Atom att2 = (Atom) temp.get(1);

							double x0 = trait.getFx();
							double y0 = trait.getFy();
							double x1 = first.getFx();
							double y1 = first.getFy();
							double x2 = last.getFx();
							double y2 = last.getFy();
							double a = y2 - y1;
							double b = x1 - x2;
							double c = y1 * x2 - x1 * y2;
							double a2 = -b;
							double b2 = a;
							double c2 = b * x0 - a * y0;
							traitAtomChaine(tempo, x0, y0, a2, b2, c2, a, b, c);

							double OAx = tempo.getFx() - trait.getFx();
							double OAy = tempo.getFy() - trait.getFy();

							double resx = trait.getFx() + Math.cos(0.69) * OAx
									- Math.sin(0.69) * OAy;
							double resy = trait.getFy() + Math.sin(0.69) * OAx
									+ Math.cos(0.69) * OAy;

							double resx2 = trait.getFx() + Math.cos(-0.69)
									* OAx - Math.sin(-0.69) * OAy;
							double resy2 = trait.getFy() + Math.sin(-0.69)
									* OAx + Math.cos(-0.69) * OAy;

							Atom ato1 = (Atom) temp2.get(0);
							Atom ato2 = (Atom) temp2.get(1);

							double angle1 = angleOld2(att, trait, att2);
							double angle2 = angleNew(resx, resy, trait, resx2,
									resy2);

							if ((angle1 < 0) && (angle2 < 0) || (angle1 > 0)
									&& (angle2 > 0)) {
								if (!att.getFlag()) {
									att.setFx(resx);
									att.setFy(resy);
									if ((att.getType() == 5)
											|| (att.getType() == 6)) {
										att.setFlag(true);
									}
								}
								if (!att2.getFlag()) {
									att2.setFx(resx2);
									att2.setFy(resy2);
									if ((att2.getType() == 5)
											|| (att2.getType() == 6)) {
										att2.setFlag(true);
									}
								}

							} else {
								if (!att.getFlag()) {
									att.setFx(resx2);
									att.setFy(resy2);
									if ((att.getType() == 5)
											|| (att.getType() == 6)) {
										att.setFlag(true);
									}
								}
								if (!att2.getFlag()) {
									att2.setFx(resx);
									att2.setFy(resy);
									if ((att2.getType() == 5)
											|| (att2.getType() == 6)) {
										att2.setFlag(true);
									}
								}

							}

						}

					}

				}

				/* traitment of the first and the last atom */

				LinkedList nt = first.getNeig();
				LinkedList nt2 = new LinkedList();

				if (nt.size() != 0) {
					for (int i = 0; i < nt.size(); i++) {
						Atom aaa = (Atom) nt.get(i);
						if (atomId.indexOf(aaa.getId()) == -1) {
							if (aaa.getType() == 6) {
								aaa.setFlag(false);
							}
							nt2.add(aaa);
						}
					}
					if (nt2.size() != 0) {
						if ((nt2.size() == 2) || (nt2.size() == 3))
							traitAtomC(first, getAtom((Integer) atomId
									.get(atomId.size() - 2)));
						if (nt2.size() == 1) {
							Atom trait = (Atom) nt2.get(0);
							intersectionImpair3(trait, first,
									getAtom((Integer) atomId
											.get(atomId.size() - 2)),
									getAtom((Integer) atomId
											.get(atomId.size() - 3)));
						}
					}
				}

				first.setFlag(true);

				Atom ult = getAtom((Integer) atomId.get(0));

				LinkedList nnt = ult.getNeig();
				LinkedList nnt2 = new LinkedList();

				if (nnt.size() != 0) {
					for (int i = 0; i < nnt.size(); i++) {
						Atom aaa = (Atom) nnt.get(i);
						if (atomId.indexOf(aaa.getId()) == -1) {
							if (aaa.getType() == 6) {
								aaa.setFlag(false);
							}
							nnt2.add(aaa);
						}
					}
					if (nnt2.size() != 0) {
						if ((nnt2.size() == 2) || (nnt2.size() == 3))
							traitAtomC(ult, getAtom((Integer) atomId.get(1)));
						if (nnt2.size() == 1) {
							Atom trait = (Atom) nnt2.get(0);
							intersectionImpair3(trait, ult,
									getAtom((Integer) atomId.get(1)),
									getAtom((Integer) atomId.get(2)));
						}
					}
				}

			}

		}

	}

	private Atom getBig(Chaine c) {
		LinkedList idAtoms = c.getAtoms();
		Atom res = getAtom((Integer) idAtoms.get(1));
		int nn = res.getNeig().size();

		for (int i = 2; i < idAtoms.size() - 1; i++) {
			Atom temp = getAtom((Integer) idAtoms.get(i));
			if (temp.getNeig().size() > nn) {
				nn = temp.getNeig().size();
				res = temp;
			}

		}

		return res;

	}

	private void traitAtomChaine(Atom trait, double xo, double yo, double a,
			double b, double c, int type) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double r = 1.2;

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

		double solx;
		double soly;

		if (type == 1) {
			if (sol1x <= sol2x) {
				solx = sol1x;
				soly = sol1y;
			} else {
				solx = sol2x;
				soly = sol2y;
			}
		} else {
			if (sol1x >= sol2x) {
				solx = sol1x;
				soly = sol1y;
			} else {
				solx = sol2x;
				soly = sol2y;
			}
		}

		trait.setFx(solx);
		trait.setFy(soly);
	}

	private void traitAtomChaine(Atom trait, double xo, double yo, double a,
			double b, double c, double a2, double b2, double c2) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double r = 1.2;

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

		double solx;
		double soly;

		double res1 = a2 * sol1x + b2 * sol1y + c2;
		double res2 = a2 * sol2x + b2 * sol2y + c2;

		if (Math.abs(res1) > Math.abs(res2)) {
			solx = sol1x;
			soly = sol1y;
		} else {
			solx = sol2x;
			soly = sol2y;

		}

		trait.setFx(solx);
		trait.setFy(soly);
	}

	public void traitAtomChaine2(Atom trait, Atom ref1, Atom ref2, Atom ref3) {

		double xo1 = ref1.getFx();
		double yo1 = ref1.getFy();
		double xo = ref2.getFx();
		double yo = ref2.getFy();

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		double r1 = 1.91;
		double r = 1.2;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo;
		double b2 = -2.0 * yo;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double a = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
		double C2 = xo * xo + yo * yo - r * r;
		double C12 = C1 - C2;
		double b, c;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);

		if (bpga) {
			b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = b * b - 4.0 * a * c;

		solx1 = -0.5 * (b + 1.0 * Math.sqrt(Math.abs(delta))) / a;
		solx2 = -0.5 * (b + -1.0 * Math.sqrt(Math.abs(delta))) / a;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			solx1 = -(b12 * solx1 + C12) / a12;
			soly2 = solx2;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		double solx;
		double soly;

		double angle1 = angle(ref1, ref2, ref3);
		double angle2 = angleNew(ref1.getFx(), ref1.getFy(), ref2, solx1, soly1);

		if (((angle1 > 0) && (angle2 < 0)) || ((angle1 < 0) && (angle2 > 0))) {
			solx = solx1;
			soly = soly1;
		} else {
			solx = solx2;
			soly = soly2;
		}

		if (!trait.getFlag()) {
			trait.setFx(solx);
			trait.setFy(soly);

		}

	}

	public double angle(Atom a1, Atom a2, Atom a3) {
		double d;

		double ux = a1.getFx() - a2.getFx();
		double uy = a1.getFy() - a2.getFy();
		double vx = a3.getFx() - a2.getFx();
		double vy = a3.getFy() - a2.getFy();

		double nu = Math.sqrt(ux * ux + uy * uy);
		double nv = Math.sqrt(vx * vx + vy * vy);

		d = Math.acos((ux * vx + uy * vy) / (nu * nv));
		if ((ux * vy - uy * vx) < 0.0)
			d = -d;

		return d;

	}

	public double angleOld(Atom a1, Atom a2, Atom a3) {
		double d;

		double ux = a1.getFx() - a2.getX();
		double uy = a1.getFy() - a2.getY();
		double vx = a3.getFx() - a2.getX();
		double vy = a3.getFy() - a2.getY();

		double nu = Math.sqrt(ux * ux + uy * uy);
		double nv = Math.sqrt(vx * vx + vy * vy);

		d = Math.acos((ux * vx + uy * vy) / (nu * nv));
		if ((ux * vy - uy * vx) < 0.0)
			d = -d;

		return d;

	}

	public double angleOld2(Atom a1, Atom a2, Atom a3) {
		double d;

		double ux = a1.getX() - a2.getX();
		double uy = a1.getY() - a2.getY();
		double vx = a3.getX() - a2.getX();
		double vy = a3.getY() - a2.getY();

		double nu = Math.sqrt(ux * ux + uy * uy);
		double nv = Math.sqrt(vx * vx + vy * vy);

		d = Math.acos((ux * vx + uy * vy) / (nu * nv));
		if ((ux * vy - uy * vx) < 0.0)
			d = -d;

		return d;

	}

	public double angleNew(double x1, double y1, Atom a2, double x2, double y2) {
		double d;

		double ux = x1 - a2.getFx();
		double uy = y1 - a2.getFy();
		double vx = x2 - a2.getFx();
		double vy = y2 - a2.getFy();

		double nu = Math.sqrt(ux * ux + uy * uy);
		double nv = Math.sqrt(vx * vx + vy * vy);

		d = Math.acos((ux * vx + uy * vy) / (nu * nv));
		if ((ux * vy - uy * vx) < 0.0)
			d = -d;

		return d;

	}

	public final double angle3(Atom a1, Atom a2, Atom a3) {

		double c1[] = new double[2];
		double c2[] = new double[2];
		double c3[] = new double[2];

		c1[0] = a1.getFx();
		c1[1] = a1.getFy();

		c2[0] = a2.getFx();
		c2[1] = a2.getFy();

		c3[0] = a3.getFx();
		c3[1] = a3.getFy();

		double[] v1 = new double[2];
		double[] v2 = new double[2];
		double angle;

		for (int i = 0; i < 2; i++) {
			v1[i] = c1[i] - c2[i];
			v2[i] = c3[i] - c2[i];
		}

		double na = Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1]);
		double nb = Math.sqrt(v2[0] * v2[0] + v2[1] * v2[1]);
		double c = (v1[0] * v2[0] + v1[1] * v2[1]) / (na * nb);
		double s = (v1[0] * v2[1]) - (v1[1] * v2[0]);
		double sign;
		if (s < 0) {
			sign = -1;
		} else
			sign = 1;

		angle = sign * Math.acos(c);

		return angle;

	}

	private void intersection(Atom first, Atom last, int n, Atom trait) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = first.getFx();
		double y1 = first.getFy();
		double x2 = last.getFx();
		double y2 = last.getFy();
		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;
		double r = n * 2.77;

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

		double solx;
		double soly;

		if (a <= 0) {
			if (sol1y < sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		} else {
			if (sol1y > sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		}

		trait.setFx(solx);
		trait.setFy(soly);

	}

	private void traitAtomRing(Atom trait, Atom ref, double x0, double y0) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x2 = ref.getFx();
		double y2 = ref.getFy();
		double a = y2 - y0;
		double b = x0 - x2;
		double c = y0 * x2 - x0 * y2;

		double xo = ref.getFx();
		double yo = ref.getFy();
		double r = 1.2;

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

		double solx;
		double soly;

		double d1 = Math.sqrt((x0 - sol1x) * (x0 - sol1x) + (y0 - sol1y)
				* (y0 - sol1y));
		double d2 = Math.sqrt((x0 - sol2x) * (x0 - sol2x) + (y0 - sol2y)
				* (y0 - sol2y));

		if (d1 < d2) {
			solx = sol2x;
			soly = sol2y;
		} else {
			solx = sol1x;
			soly = sol1y;
		}

		trait.setFx(solx);
		trait.setFy(soly);

	}

	public boolean isInsideId(Atom at, LinkedList l) {

		for (int i = 0; i < l.size(); i++) {
			int id = (Integer) l.get(i);
			if (at.getId() == id)
				return true;
		}
		return false;
	}

	public void traitRing(Atom dep, Ring ring) {
		boolean conj = false;
		if (ring.getFlag())
			return;
		if (ring.getConjugate().size() != 0) {

			LinkedList atoms = ring.getAtoms();
			for (int i = 0; i < atoms.size(); i++) {
				Atom att = getAtom((Integer) atoms.get(i));
				if (att.getFlag()) {
					dep = att;
					conj = true;
				}
			}

		}
		boolean conj2 = false;
		if (ring.getConjugate().size() != 0) {
			conj2 = true;
		}

		int n = ring.getLenght();
		ring.setFlag(true);
		Atom[] ringAtoms = new Atom[n];
		LinkedList temp = ring.getAtoms();

		int pos = 0;

		for (int i = 0; i < temp.size(); i++) {
			if ((Integer) temp.get(i) == dep.getId()) {
				pos = i;
				break;
			}
		}

		ringAtoms[0] = getAtom((Integer) temp.get(pos));

		Atom a1;
		Atom a2;
		if (pos == 0) {
			a1 = getAtom((Integer) temp.get(temp.size() - 1));

			a2 = getAtom((Integer) temp.get(pos + 1));
		} else {
			if (pos == temp.size() - 1) {
				a1 = getAtom((Integer) temp.get(pos - 2));
				a2 = getAtom((Integer) temp.get(0));
			} else {
				a1 = getAtom((Integer) temp.get(pos - 1));
				a2 = getAtom((Integer) temp.get(pos + 1));
			}

		}

		int sens = 0;

		double angle = angleOld2(a1, dep, a2);

		if (angle >= 0)
			sens = 1;
		else
			sens = -1;

		int cont = 1;

		if (sens == 1) {
			for (int i = pos + 1; i < temp.size(); i++) {
				ringAtoms[cont] = getAtom((Integer) temp.get(i));
				cont++;
			}
			for (int i = 0; i < pos; i++) {
				ringAtoms[cont] = getAtom((Integer) temp.get(i));
				cont++;
			}

		}

		if (sens == -1) {
			if (pos != 0) {
				for (int i = pos - 1; i >= 0; i--) {
					ringAtoms[cont] = getAtom((Integer) temp.get(i));
					cont++;
				}
			}
			for (int i = temp.size() - 1; i > pos; i--) {
				ringAtoms[cont] = getAtom((Integer) temp.get(i));
				cont++;
			}

		}

		if (n % 2 == 0) {

			double r = 0.6 / (Math.sin(Math.toRadians(180 / n)));
			double xo, yo;

			Atom prede2 = getPredecesor(ring);

			if (!conj) {
				if (prede2 != null) {

					Atom prede = prede2.clone();
					double nx = 2.0 * ringAtoms[0].getFx() - prede2.getFx();
					double ny = 2.0 * ringAtoms[0].getFy() - prede2.getFy();
					prede.setFx(nx);
					prede.setFy(ny);

					xo = getCentreX(ringAtoms[0], prede, r);
					yo = getCentreY(ringAtoms[0], prede, r);

				} else {
					xo = getCentreX(ringAtoms[0], ringAtoms[n / 2], r);
					yo = getCentreY(ringAtoms[0], ringAtoms[n / 2], r);

				}

			} else {

				LinkedList rings = ring.getConjugate();
				Ring conjugate = null;
				for (int i = 0; i < rings.size(); i++) {
					Ring rr = (Ring) rings.get(i);
					if (rr.getFlag())
						conjugate = rr;
					continue;
				}

				Atom ref1 = null;
				Atom ref2 = null;

				LinkedList tt = new LinkedList();

				int cc = 0;

				for (int i = 0; i < ringAtoms.length - 1; i++) {
					if ((ringAtoms[i].getFlag())
							&& (ringAtoms[i + 1].getFlag())) {
						tt.add(ringAtoms[i]);
						cc++;
						tt.add(ringAtoms[i + 1]);
						cc++;
					}
					if ((ringAtoms[0].getFlag())
							&& (ringAtoms[ringAtoms.length - 1].getFlag())) {
						tt.add(ringAtoms[i]);
						cc++;
						tt.add(ringAtoms[ringAtoms.length - 1]);
						cc++;
					}
					if (cc == 2)
						break;
				}

				ref1 = (Atom) tt.get(0);
				ref2 = (Atom) tt.get(1);

				xo = getCenterConjX(ref1, ref2, r, conjugate);
				yo = getCenterConjY(ref1, ref2, r, conjugate);

			}

			if (ring.isAromatic()) {
				ring.setCenterX(xo);

				ring.setCenterY(yo);
			}

			intersectionFirst(dep, xo, yo, r, ringAtoms, sens);

			for (int i = 1; i < n - 1; i++) {
				Atom atom = ringAtoms[i];

				intersectionRing(atom, xo, yo, r, ringAtoms, i);

				atom.setFlag(true);
			}
			for (int i = 0; i < n; i++) {
				Atom atom = ringAtoms[i];

				LinkedList neig = atom.getNeig();
				LinkedList temp1 = new LinkedList();
				LinkedList temp2 = new LinkedList();
				for (int s = 0; s < neig.size(); s++) {
					Atom aa = (Atom) neig.get(s);
					if (temp.indexOf(aa.getId()) == -1) {
						temp1.add(aa);
					} else
						temp2.add(aa);
				}
				if (temp1.size() == 1) {
					Atom att = (Atom) temp1.get(0);
					if ((!att.getFlag()) || (att.getType() == 6)) {

						traitAtomRing(att, atom, xo, yo);

						if (!conj2) {
							traitAtom(att);
						} else {
							if (!isConjugate(att, ring)) {
								traitAtom(att);
							}
						}
						if ((att.getType() == 5) || (att.getType() == 6)) {
							att.setFlag(true);
						}

					}
				}

				if (temp1.size() == 2) {

					Atom att = (Atom) temp1.get(0);
					Atom att2 = (Atom) temp1.get(1);

					Atom tempo = att.clone();
					traitAtomRing(tempo, atom, xo, yo);

					double OAx = tempo.getFx() - atom.getFx();
					double OAy = tempo.getFy() - atom.getFy();

					double resx = atom.getFx() + Math.cos(0.69) * OAx
							- Math.sin(0.69) * OAy;
					double resy = atom.getFy() + Math.sin(0.69) * OAx
							+ Math.cos(0.69) * OAy;

					double resx2 = atom.getFx() + Math.cos(-0.69) * OAx
							- Math.sin(-0.69) * OAy;
					double resy2 = atom.getFy() + Math.sin(-0.69) * OAx
							+ Math.cos(-0.69) * OAy;

					double d1 = Math.sqrt((resx - att.getFx())
							* (resx - att.getFx()) + (resy - att.getFy())
							* (resy - att.getFy()));
					double d2 = Math.sqrt((resx2 - att.getFx())
							* (resx2 - att.getFx()) + (resy2 - att.getFy())
							* (resy2 - att.getFy()));

					if (d1 > d2) {
						if ((!att.getFlag()) || (att.getType() == 6)) {
							att.setFx(resx);
							att.setFy(resy);

							if (!conj2) {
								traitAtom(att);
							} else {
								if (!isConjugate(att, ring)) {
									traitAtom(att);
								}
							}
							if ((att.getType() == 5) || (att.getType() == 6)) {
								att.setFlag(true);
							}
						}
						if ((!att2.getFlag()) || (att2.getType() == 6)) {

							att2.setFx(resx2);
							att2.setFy(resy2);
							if (!conj2) {
								traitAtom(att2);
							} else {
								if (!isConjugate(att2, ring)) {
									traitAtom(att2);
								}
							}
							if ((att2.getType() == 5) || (att2.getType() == 6)) {
								att2.setFlag(true);
							}
						}
					} else {
						if ((!att.getFlag()) || (att.getType() == 6)) {

							att.setFx(resx2);
							att.setFy(resy2);
							if (!conj2) {
								traitAtom(att);
							} else {
								if (!isConjugate(att, ring)) {
									traitAtom(att);
								}
							}
							if ((att.getType() == 5) || (att.getType() == 6)) {
								att.setFlag(true);
							}
						}
						if ((!att2.getFlag()) || (att2.getType() == 6)) {

							att2.setFx(resx);
							att2.setFy(resy);
							if (!conj2) {
								traitAtom(att2);
							} else {
								if (!isConjugate(att2, ring)) {
									traitAtom(att2);
								}
							}
							if ((att2.getType() == 5) || (att2.getType() == 6)) {
								att2.setFlag(true);
							}
						}
					}

				}

			}

		} else {

			double r = 0.6 / (Math.sin(Math.toRadians(180 / n)));
			double xo;
			double yo;
			if (!conj) {
				Atom prede2 = getPredecesor(ring);

				if (prede2 != null) {

					Atom prede = prede2.clone();
					double nx = 2.0 * ringAtoms[0].getFx() - prede2.getFx();
					double ny = 2.0 * ringAtoms[0].getFy() - prede2.getFy();
					prede.setFx(nx);
					prede.setFy(ny);

					xo = getCentreX(ringAtoms[0], prede, r);
					yo = getCentreY(ringAtoms[0], prede, r);

				} else {

					xo = getCentre2X(ringAtoms, r);
					yo = getCentre2Y(ringAtoms, r);
				}
			} else {

				LinkedList rings = ring.getConjugate();
				Ring conjugate = null;
				for (int i = 0; i < rings.size(); i++) {
					Ring rr = (Ring) rings.get(i);
					if (rr.getFlag())
						conjugate = rr;
					continue;
				}

				Atom ref1 = null;
				Atom ref2 = null;

				LinkedList tt = new LinkedList();

				int cc = 0;

				for (int i = 0; i < ringAtoms.length; i++) {
					if (ringAtoms[i].getFlag()) {
						tt.add(ringAtoms[i]);
						cc++;
					}
					if (cc == 2)
						continue;
				}

				ref1 = (Atom) tt.get(0);
				ref2 = (Atom) tt.get(1);

				xo = getCenterConjX(ref1, ref2, r, conjugate);
				yo = getCenterConjY(ref1, ref2, r, conjugate);

			}

			intersectionFirst(dep, xo, yo, r, ringAtoms, sens);

			for (int i = 1; i < n - 1; i++) {
				Atom atom = ringAtoms[i];

				intersectionRing(atom, xo, yo, r, ringAtoms, i);

			}

			for (int i = 0; i < n; i++) {
				Atom atom = ringAtoms[i];

				LinkedList neig = atom.getNeig();
				LinkedList temp1 = new LinkedList();
				LinkedList temp2 = new LinkedList();
				for (int s = 0; s < neig.size(); s++) {
					Atom aa = (Atom) neig.get(s);
					if (temp.indexOf(aa.getId()) == -1) {
						temp1.add(aa);
					} else
						temp2.add(aa);
				}
				if (temp1.size() == 1) {
					Atom att = (Atom) temp1.get(0);

					if ((!att.getFlag()) || (att.getType() == 6)) {

						traitAtomRing(att, atom, xo, yo);
						if (!conj2) {
							traitAtom(att);
						} else {
							if (!isConjugate(att, ring)) {
								traitAtom(att);
							}
						}
						if ((att.getType() == 5) || (att.getType() == 6)) {
							att.setFlag(true);
						}
					}

				}

				if (temp1.size() == 2) {
					Atom att = (Atom) temp1.get(0);
					Atom att2 = (Atom) temp1.get(1);

					Atom tempo = att.clone();
					traitAtomRing(tempo, atom, xo, yo);

					double OAx = tempo.getFx() - atom.getFx();
					double OAy = tempo.getFy() - atom.getFy();

					double resx = atom.getFx() + Math.cos(0.69) * OAx
							- Math.sin(0.69) * OAy;
					double resy = atom.getFy() + Math.sin(0.69) * OAx
							+ Math.cos(0.69) * OAy;

					double resx2 = atom.getFx() + Math.cos(-0.69) * OAx
							- Math.sin(-0.69) * OAy;
					double resy2 = atom.getFy() + Math.sin(-0.69) * OAx
							+ Math.cos(-0.69) * OAy;

					double d1 = Math.sqrt((resx - att.getFx())
							* (resx - att.getFx()) + (resy - att.getFy())
							* (resy - att.getFy()));
					double d2 = Math.sqrt((resx2 - att.getFx())
							* (resx2 - att.getFx()) + (resy2 - att.getFy())
							* (resy2 - att.getFy()));

					if (d1 > d2) {
						if ((!att.getFlag()) || (att.getType() == 6)) {
							att.setFx(resx);
							att.setFy(resy);

							if (!conj2) {
								traitAtom(att);
							} else {
								if (!isConjugate(att, ring)) {
									traitAtom(att);
								}
							}
							if ((att.getType() == 5) || (att.getType() == 6)) {
								att.setFlag(true);
							}
						}
						if ((!att2.getFlag()) || (att2.getType() == 6)) {

							att2.setFx(resx2);
							att2.setFy(resy2);
							if (!conj2) {
								traitAtom(att2);
							} else {
								if (!isConjugate(att2, ring)) {
									traitAtom(att2);
								}
							}
							if ((att2.getType() == 5) || (att2.getType() == 6)) {
								att2.setFlag(true);
							}
						}
					} else {
						if ((!att.getFlag()) || (att.getType() == 6)) {

							att.setFx(resx2);
							att.setFy(resy2);
							if (!conj2) {
								traitAtom(att);
							} else {
								if (!isConjugate(att, ring)) {
									traitAtom(att);
								}
							}
							if ((att.getType() == 5) || (att.getType() == 6)) {
								att.setFlag(true);
							}
						}
						if ((!att2.getFlag()) || (att2.getType() == 6)) {

							att2.setFx(resx);
							att2.setFy(resy);
							if (!conj2) {
								traitAtom(att2);
							} else {
								if (!isConjugate(att2, ring)) {
									traitAtom(att2);
								}
							}
							if ((att2.getType() == 5) || (att2.getType() == 6)) {
								att2.setFlag(true);
							}
						}
					}

				}

			}

		}

		double[] coord = new double[4];
		coord[0] = 9.9e19;
		coord[1] = -9.9e19;
		coord[2] = 9.9e19;
		coord[3] = -9.9e19;

		for (int i = 0; i < ringAtoms.length; i++) {
			if (ringAtoms[i].getFx() < coord[0]) {
				coord[0] = ringAtoms[i].getFx();
			}
			if (ringAtoms[i].getFx() > coord[1]) {
				coord[1] = ringAtoms[i].getFx();
			}
			if (ringAtoms[i].getFy() < coord[2]) {
				coord[2] = ringAtoms[i].getFy();
			}
			if (ringAtoms[i].getFy() > coord[3]) {
				coord[3] = ringAtoms[i].getFy();
			}
		}

		ring.setCoord(coord);

	}

	public boolean isConjugate(Atom atom, Ring r) {

		LinkedList l = r.getConjugate();
		for (int i = 0; i < l.size(); i++) {
			Ring rr = (Ring) l.get(i);
			LinkedList aa = rr.getAtoms();

			for (int j = 0; j < aa.size(); j++) {
				if (atom.getId() == (Integer) aa.get(j))
					return true;
			}

		}
		return false;
	}

	public void intersectionFirst(Atom dep, double xo, double yo, double r,
			Atom[] tab, int sens) {

		double xo1 = dep.getFx();
		double yo1 = dep.getFy();

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		double r1 = 1.2;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo;
		double b2 = -2.0 * yo;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double a = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
		double C2 = xo * xo + yo * yo - r * r;
		double C12 = C1 - C2;
		double b, c;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);

		if (bpga) {
			b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = b * b - 4.0 * a * c;

		solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
		solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			solx1 = -(b12 * solx1 + C12) / a12;
			soly2 = solx2;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		Atom at1 = tab[1];
		Atom at2 = tab[tab.length - 1];

		if ((!at1.getFlag()) && (!at2.getFlag())) {
			double angle1 = angleNew(solx1, soly1, dep, solx2, soly2);
			double angle2 = angleOld2(at1, dep, at2);

			if (angle1 < 0) {
				if (!at1.getFlag()) {
					at1.setFx(solx1);
					at1.setFy(soly1);
					at1.setFlag(true);
				}
				if (!at2.getFlag()) {
					at2.setFx(solx2);
					at2.setFy(soly2);
					at2.setFlag(true);
				}
				dep.setFlag(true);
			} else {
				if (!at1.getFlag()) {
					at1.setFx(solx2);
					at1.setFy(soly2);
					at1.setFlag(true);
				}
				if (!at2.getFlag()) {
					at2.setFx(solx1);
					at2.setFy(soly1);
					at2.setFlag(true);
				}
				dep.setFlag(true);
			}
		}

		else {

			Atom fixe = null;
			Atom Nfixe = null;
			if (at1.getFlag()) {
				fixe = at1;
				Nfixe = at2;
			} else {
				fixe = at2;
				Nfixe = at1;
			}

			double dd1 = Math.sqrt((fixe.getFx() - solx1)
					* (fixe.getFx() - solx1) + (fixe.getFy() - soly1)
					* (fixe.getFy() - soly1));
			double dd2 = Math.sqrt((fixe.getFx() - solx2)
					* (fixe.getFx() - solx2) + (fixe.getFy() - soly2)
					* (fixe.getFy() - soly2));

			if (!Nfixe.getFlag()) {
				if (dd1 < dd2) {
					Nfixe.setFx(solx2);
					Nfixe.setFy(soly2);
					Nfixe.setFlag(true);
				} else {
					Nfixe.setFx(solx1);
					Nfixe.setFy(soly1);
					Nfixe.setFlag(true);
				}

			}

		}

	}

	public void intersectionRing(Atom dep, double xo, double yo, double r,
			Atom[] tab, int pos) {

		double xo1 = dep.getFx();
		double yo1 = dep.getFy();

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		double r1 = 1.2;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo;
		double b2 = -2.0 * yo;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double a = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r1 * r1;
		double C2 = xo * xo + yo * yo - r * r;
		double C12 = C1 - C2;
		double b, c;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);

		if (bpga) {
			b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = b * b - (4.0 * a * c);

		solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
		solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			solx1 = -(b12 * solx1 + C12) / a12;
			soly2 = solx2;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		double solx;
		double soly;

		Atom at1 = tab[pos - 1];
		Atom at2 = tab[pos + 1];

		if (Double.isNaN(solx1)) {

			Atom con = dep.clone();

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x2 = dep.getFx();
			double y2 = dep.getFy();
			double aaa = y2 - yo;
			double bbb = xo - x2;
			double ccc = yo * x2 - xo * y2;

			double rrr = 1.2;

			if (b != 0.0) {
				double u = aaa * (ccc + bbb * yo) - bbb * bbb * xo;
				double v = aaa * aaa + bbb * bbb;
				double w = ccc + bbb * yo;
				double deltap = u * u - v
						* (bbb * bbb * (xo * xo - rrr * rrr) + w * w);
				if (deltap >= 0.0) {
					sol1x = (-u + 1.0 * Math.sqrt(deltap)) / v;
					sol1y = -(aaa * sol1x + ccc) / bbb;
					sol2x = (-u + -1.0 * Math.sqrt(deltap)) / v;
					sol2y = -(aaa * sol2x + ccc) / bbb;
				}

			} else {
				sol1x = -ccc / aaa;
				double u = sol1x - xo;
				double v = rrr * rrr - u * u;
				if (v >= 0.0) {
					sol1y = yo + 1.0 * Math.sqrt(v);
					sol2y = yo + -1.0 * Math.sqrt(v);
					sol2x = sol1x;
				}

			}

			double d1 = Math.sqrt((dep.getFx() - sol1x) * (dep.getFx() - sol1x)
					+ (dep.getFy() - sol1y) * (dep.getFy() - sol1y));
			double d2 = Math.sqrt((dep.getFx() - sol2x) * (dep.getFx() - sol2x)
					+ (dep.getFy() - sol2y) * (dep.getFy() - sol2y));

			if (d1 < d2) {
				con.setFx(sol1x);
				con.setFy(sol1y);
			} else {
				con.setFx(sol2x);
				con.setFy(sol2y);
			}

			double xoo1 = con.getFx();
			double yoo1 = con.getFy();

			double rr1 = 1.2;

			double aa1 = -2.0 * xoo1;
			double bb1 = -2.0 * yoo1;
			double aa2 = -2.0 * xo;
			double bb2 = -2.0 * yo;
			double aa12 = aa1 - aa2;
			double bb12 = bb1 - bb2;
			double aa = aa12 * aa12 + bb12 * bb12;

			double Cc1 = xoo1 * xoo1 + yoo1 * yoo1 - rr1 * rr1;
			double Cc2 = xo * xo + yo * yo - r * r;
			double Cc12 = Cc1 - Cc2;
			double bb, cc;
			boolean bpga2 = Math.abs(bb12) >= Math.abs(aa12);

			if (bpga2) {
				bb = 2.0 * aa12 * Cc12 + bb12 * (bb1 * aa2 - aa1 * bb2);
				cc = Cc12 * Cc12 + bb12 * (bb1 * Cc2 - bb2 * Cc1);
			} else {
				bb = 2.0 * bb12 * Cc12 + aa12 * (aa1 * bb2 - bb1 * aa2);
				cc = Cc12 * Cc12 + aa12 * (aa1 * Cc2 - aa2 * Cc1);
			}
			double delta2 = bb * bb - (4.0 * aa * cc);

			solx1 = -0.5 * (bb + 1.0 * Math.sqrt(delta2)) / aa;
			solx2 = -0.5 * (bb + -1.0 * Math.sqrt(delta2)) / aa;

			if (bpga2) {
				soly1 = -(aa12 * solx1 + Cc12) / bb12;
				soly2 = -(aa12 * solx2 + Cc12) / bb12;
			} else {
				soly1 = solx1;
				solx1 = -(bb12 * solx1 + Cc12) / aa12;
				soly2 = solx2;
				solx2 = -(bb12 * solx2 + Cc12) / aa12;
			}

		}

		double dx1 = Math.abs(solx1 - at1.getFx());
		double dx2 = Math.abs(solx2 - at1.getFx());
		double dy1 = Math.abs(soly1 - at1.getFy());
		double dy2 = Math.abs(soly2 - at1.getFy());

		if ((dx1 < dx2) && (dy1 < dy2)) {
			solx = solx2;
			soly = soly2;
		} else {
			solx = solx1;
			soly = soly1;
		}

		if (!at2.getFlag()) {

			at2.setFx(solx);
			at2.setFy(soly);
			at2.setFlag(true);
		}

	}

	private double getCentreX(Atom first, Atom last, double r) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = first.getFx();
		double y1 = first.getFy();
		double x2 = last.getFx();
		double y2 = last.getFy();
		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;

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

		double solx;
		double soly;

		if (a <= 0) {
			if (sol1y < sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		} else {
			if (sol1y > sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		}

		return solx;

	}

	private double getCenterConjX(Atom ref1, Atom ref2, double r, Ring ring) {

		double xo1 = ref1.getFx();
		double yo1 = ref1.getFy();
		double xo2 = ref2.getFx();
		double yo2 = ref2.getFy();

		double solx = 0;
		double soly = 0;

		if (Math.sqrt((xo1 - xo2) * (xo1 - xo2) + (yo1 - yo2) * (yo1 - yo2)) - 1.2 < 0.1) {

			double a1 = -2.0 * xo1;
			double b1 = -2.0 * yo1;
			double a2 = -2.0 * xo2;
			double b2 = -2.0 * yo2;
			double a12 = a1 - a2;
			double b12 = b1 - b2;
			double a = a12 * a12 + b12 * b12;

			double C1 = xo1 * xo1 + yo1 * yo1 - r * r;
			double C2 = xo2 * xo2 + yo2 * yo2 - r * r;
			double C12 = C1 - C2;
			double b, c;
			boolean bpga = Math.abs(b12) >= Math.abs(a12);
			if (bpga) {
				b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
				c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
			} else {
				b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
				c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
			}
			double delta = b * b - 4.0 * a * c;

			double solx1 = 0;
			double soly1 = 0;
			double solx2 = 0;
			double soly2 = 0;

			solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
			solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

			if (bpga) {
				soly1 = -(a12 * solx1 + C12) / b12;
				soly2 = -(a12 * solx2 + C12) / b12;
			} else {
				soly1 = solx1;
				soly2 = solx2;
				solx1 = -(b12 * solx1 + C12) / a12;
				solx2 = -(b12 * solx2 + C12) / a12;
			}

			if (ringContains(ring, solx1, soly1) == 1) {
				solx = solx2;
				soly = soly2;
			}

			else {
				solx = solx1;
				soly = soly1;
			}

		} else {

			Atom con = ref2.clone();

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x2 = ref2.getFx();
			double y2 = ref2.getFy();
			double x1 = ref1.getFx();
			double y1 = ref1.getFy();
			double aaa = y2 - y1;
			double bbb = x1 - x2;
			double ccc = y1 * x2 - x1 * y2;

			double rrr = 1.2;

			if (bbb != 0.0) {
				double u = aaa * (ccc + bbb * y1) - bbb * bbb * x1;
				double v = aaa * aaa + bbb * bbb;
				double w = ccc + bbb * y1;
				double deltap = u * u - v
						* (bbb * bbb * (x1 * x1 - rrr * rrr) + w * w);
				if (deltap >= 0.0) {
					sol1x = (-u + 1.0 * Math.sqrt(deltap)) / v;
					sol1y = -(aaa * sol1x + ccc) / bbb;
					sol2x = (-u + -1.0 * Math.sqrt(deltap)) / v;
					sol2y = -(aaa * sol2x + ccc) / bbb;
				}

			} else {
				sol1x = -ccc / aaa;
				double u = sol1x - x1;
				double v = rrr * rrr - u * u;
				if (v >= 0.0) {
					sol1y = y1 + 1.0 * Math.sqrt(v);
					sol2y = y1 + -1.0 * Math.sqrt(v);
					sol2x = sol1x;
				}

			}

			double d1 = Math.sqrt((ref2.getFx() - sol1x)
					* (ref2.getFx() - sol1x) + (ref2.getFy() - sol1y)
					* (ref2.getFy() - sol1y));
			double d2 = Math.sqrt((ref2.getFx() - sol2x)
					* (ref2.getFx() - sol2x) + (ref2.getFy() - sol2y)
					* (ref2.getFy() - sol2y));

			if (d1 < d2) {
				con.setFx(sol1x);
				con.setFy(sol1y);
			} else {
				con.setFx(sol2x);
				con.setFy(sol2y);
			}

			xo2 = con.getFx();
			yo2 = con.getFy();

			double a1 = -2.0 * xo1;
			double b1 = -2.0 * yo1;
			double a2 = -2.0 * xo2;
			double b2 = -2.0 * yo2;
			double a12 = a1 - a2;
			double b12 = b1 - b2;
			double a = a12 * a12 + b12 * b12;

			double C1 = xo1 * xo1 + yo1 * yo1 - r * r;
			double C2 = xo2 * xo2 + yo2 * yo2 - r * r;
			double C12 = C1 - C2;
			double b, c;
			boolean bpga = Math.abs(b12) >= Math.abs(a12);
			if (bpga) {
				b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
				c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
			} else {
				b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
				c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
			}
			double delta = b * b - 4.0 * a * c;

			double solx1 = 0;
			double soly1 = 0;
			double solx2 = 0;
			double soly2 = 0;

			solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
			solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

			if (bpga) {
				soly1 = -(a12 * solx1 + C12) / b12;
				soly2 = -(a12 * solx2 + C12) / b12;
			} else {
				soly1 = solx1;
				soly2 = solx2;
				solx1 = -(b12 * solx1 + C12) / a12;
				solx2 = -(b12 * solx2 + C12) / a12;
			}

			if (ringContains(ring, solx1, soly1) == 1) {
				solx = solx2;
				soly = soly2;
			}

			else {
				solx = solx1;
				soly = soly1;
			}

		}

		return solx;

	}

	public int ringContains(Ring r, double x, double y) {

		LinkedList l = r.getAtoms();
		Atom[] points = new Atom[l.size()];

		for (int i = 0; i < l.size(); i++) {
			points[i] = getAtom((Integer) l.get(i));
		}

		int i, j, c = 0;
		for (i = 0, j = points.length - 1; i < points.length; j = i++) {
			if ((((points[i].getFy() <= y) && (y < points[j].getFy())) || ((points[j]
					.getFy() <= y) && (y < points[i].getFy())))
					&& (x < (points[j].getFx() - points[i].getFx())
							* (y - points[i].getFy())
							/ (points[j].getFy() - points[i].getFy())
							+ points[i].getFx()))
				if (c == 0)
					c = 1;
				else
					c = 0;
		}
		return c;
	}

	public boolean ringContains2(Ring r, double x, double y) {

		LinkedList l = r.getAtoms();
		Atom[] points = new Atom[l.size()];

		for (int i = 0; i < l.size(); i++) {
			points[i] = getAtom((Integer) l.get(i));
		}

		int crossings = 0;
		for (int i = 0; i < points.length - 1; i++) {
			double slope = (points[i + 1].getFx() - points[i].getFx())
					/ (points[i + 1].getFy() - points[i].getFy());
			boolean cond1 = (points[i].getFy() <= y)
					&& (y < points[i + 1].getFy());
			boolean cond2 = (points[i + 1].getFy() <= y)
					&& (y < points[i].getFy());
			boolean cond3 = x < slope * (y - points[i].getFy())
					+ points[i].getFx();
			if ((cond1 || cond2) && cond3)
				crossings++;
		}
		return (crossings % 2 != 0);
	}

	private double getCenterConjY(Atom ref1, Atom ref2, double r, Ring ring) {

		double xo1 = ref1.getFx();
		double yo1 = ref1.getFy();
		double xo2 = ref2.getFx();
		double yo2 = ref2.getFy();

		double solx = 0;
		double soly = 0;

		if (Math.sqrt((xo1 - xo2) * (xo1 - xo2) + (yo1 - yo2) * (yo1 - yo2)) - 1.2 < 0.1) {

			double a1 = -2.0 * xo1;
			double b1 = -2.0 * yo1;
			double a2 = -2.0 * xo2;
			double b2 = -2.0 * yo2;
			double a12 = a1 - a2;
			double b12 = b1 - b2;
			double a = a12 * a12 + b12 * b12;

			double C1 = xo1 * xo1 + yo1 * yo1 - r * r;
			double C2 = xo2 * xo2 + yo2 * yo2 - r * r;
			double C12 = C1 - C2;
			double b, c;
			boolean bpga = Math.abs(b12) >= Math.abs(a12);
			if (bpga) {
				b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
				c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
			} else {
				b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
				c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
			}
			double delta = b * b - 4.0 * a * c;

			double solx1 = 0;
			double soly1 = 0;
			double solx2 = 0;
			double soly2 = 0;

			solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
			solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

			if (bpga) {
				soly1 = -(a12 * solx1 + C12) / b12;
				soly2 = -(a12 * solx2 + C12) / b12;
			} else {
				soly1 = solx1;
				soly2 = solx2;
				solx1 = -(b12 * solx1 + C12) / a12;
				solx2 = -(b12 * solx2 + C12) / a12;
			}

			if (ringContains(ring, solx1, soly1) == 1) {
				solx = solx2;
				soly = soly2;
			}

			else {
				solx = solx1;
				soly = soly1;
			}

		} else {

			Atom con = ref2.clone();

			double sol1x = 0;
			double sol1y = 0;
			double sol2x = 0;
			double sol2y = 0;

			double x2 = ref2.getFx();
			double y2 = ref2.getFy();
			double x1 = ref1.getFx();
			double y1 = ref1.getFy();
			double aaa = y2 - y1;
			double bbb = x1 - x2;
			double ccc = y1 * x2 - x1 * y2;

			double rrr = 1.2;

			if (bbb != 0.0) {
				double u = aaa * (ccc + bbb * y1) - bbb * bbb * x1;
				double v = aaa * aaa + bbb * bbb;
				double w = ccc + bbb * y1;
				double deltap = u * u - v
						* (bbb * bbb * (x1 * x1 - rrr * rrr) + w * w);
				if (deltap >= 0.0) {
					sol1x = (-u + 1.0 * Math.sqrt(deltap)) / v;
					sol1y = -(aaa * sol1x + ccc) / bbb;
					sol2x = (-u + -1.0 * Math.sqrt(deltap)) / v;
					sol2y = -(aaa * sol2x + ccc) / bbb;
				}

			} else {
				sol1x = -ccc / aaa;
				double u = sol1x - x1;
				double v = rrr * rrr - u * u;
				if (v >= 0.0) {
					sol1y = y1 + 1.0 * Math.sqrt(v);
					sol2y = y1 + -1.0 * Math.sqrt(v);
					sol2x = sol1x;
				}

			}

			double d1 = Math.sqrt((ref2.getFx() - sol1x)
					* (ref2.getFx() - sol1x) + (ref2.getFy() - sol1y)
					* (ref2.getFy() - sol1y));
			double d2 = Math.sqrt((ref2.getFx() - sol2x)
					* (ref2.getFx() - sol2x) + (ref2.getFy() - sol2y)
					* (ref2.getFy() - sol2y));

			if (d1 < d2) {
				con.setFx(sol1x);
				con.setFy(sol1y);
			} else {
				con.setFx(sol2x);
				con.setFy(sol2y);
			}

			xo2 = con.getFx();
			yo2 = con.getFy();

			double a1 = -2.0 * xo1;
			double b1 = -2.0 * yo1;
			double a2 = -2.0 * xo2;
			double b2 = -2.0 * yo2;
			double a12 = a1 - a2;
			double b12 = b1 - b2;
			double a = a12 * a12 + b12 * b12;

			double C1 = xo1 * xo1 + yo1 * yo1 - r * r;
			double C2 = xo2 * xo2 + yo2 * yo2 - r * r;
			double C12 = C1 - C2;
			double b, c;
			boolean bpga = Math.abs(b12) >= Math.abs(a12);
			if (bpga) {
				b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
				c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
			} else {
				b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
				c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
			}
			double delta = b * b - 4.0 * a * c;

			double solx1 = 0;
			double soly1 = 0;
			double solx2 = 0;
			double soly2 = 0;

			solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
			solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

			if (bpga) {
				soly1 = -(a12 * solx1 + C12) / b12;
				soly2 = -(a12 * solx2 + C12) / b12;
			} else {
				soly1 = solx1;
				soly2 = solx2;
				solx1 = -(b12 * solx1 + C12) / a12;
				solx2 = -(b12 * solx2 + C12) / a12;
			}

			if (ringContains(ring, solx1, soly1) == 1) {
				solx = solx2;
				soly = soly2;
			}

			else {
				solx = solx1;
				soly = soly1;
			}

		}

		return soly;

	}

	private double getCentre2X(Atom[] tab, double r) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1, y1, x2, y2;

		x1 = tab[0].getFx();
		y1 = tab[0].getFy();

		Atom prede2 = getPredecesor(tab[0]);

		Atom prede;

		if (prede2 != null) {
			prede = prede2.clone();
			double nx = 2.0 * tab[0].getFx() - prede2.getFx();
			double ny = 2.0 * tab[0].getFy() - prede2.getFy();
			prede.setFx(nx);
			prede.setFy(ny);
			x2 = prede.getFx();
			y2 = prede.getFy();

		} else {

			int pos = (tab.length + 1) / 2;

			Atom a1 = tab[pos];
			Atom a2;
			a2 = tab[pos - 1];

			x2 = (a1.getFx() + a2.getFx()) / 2;
			y2 = (a1.getFy() + a2.getFy()) / 2;
		}

		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;

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

		double solx;
		double soly;

		if (a <= 0) {
			if (sol1y < sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		} else {
			if (sol1y > sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		}

		return solx;

	}

	private double getCentre2ConjX(Atom[] tab, Atom ref1, Atom ref2, double r) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1, y1, x2, y2;

		x1 = tab[0].getFx();
		y1 = tab[0].getFy();

		Atom prede2 = getPredecesor(tab[0]);

		Atom prede;

		if (prede2 != null) {
			prede = prede2.clone();
			double nx = 2.0 * tab[0].getFx() - prede2.getFx();
			double ny = 2.0 * tab[0].getFy() - prede2.getFy();
			prede.setFx(nx);
			prede.setFy(ny);
			x2 = prede.getFx();
			y2 = prede.getFy();

		} else {

			int pos = (tab.length + 1) / 2;

			Atom a1 = tab[pos];
			Atom a2;
			a2 = tab[pos - 1];

			x2 = (a1.getFx() + a2.getFx()) / 2;
			y2 = (a1.getFy() + a2.getFy()) / 2;
		}

		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;

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

		double solx;
		double soly;

		if (a <= 0) {
			if (sol1y < sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		} else {
			if (sol1y > sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		}

		return solx;

	}

	private double getCentre2Y(Atom[] tab, double r) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = tab[0].getFx();
		double y1 = tab[0].getFy();

		double x2, y2;

		Atom prede2 = getPredecesor(tab[0]);

		Atom prede;

		if (prede2 != null) {
			prede = prede2.clone();
			double nx = 2.0 * tab[0].getFx() - prede2.getFx();
			double ny = 2.0 * tab[0].getFy() - prede2.getFy();
			prede.setFx(nx);
			prede.setFy(ny);
			x2 = prede.getFx();
			y2 = prede.getFy();

		} else {
			int pos = (tab.length + 1) / 2;

			Atom a1 = tab[pos];
			Atom a2;

			a2 = tab[pos - 1];

			x2 = (a1.getFx() + a2.getFx()) / 2;
			y2 = (a1.getFy() + a2.getFy()) / 2;
		}

		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;

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

		double solx;
		double soly;

		if (a <= 0) {
			if (sol1y < sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		} else {
			if (sol1y > sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		}

		return soly;

	}

	private double getCentreY(Atom first, Atom last, double r) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double x1 = first.getFx();
		double y1 = first.getFy();
		double x2 = last.getFx();
		double y2 = last.getFy();
		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;

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

		double solx;
		double soly;

		if (a <= 0) {
			if (sol1y < sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		} else {
			if (sol1y > sol2y) {
				soly = sol1y;
				solx = sol1x;
			} else {
				soly = sol2y;
				solx = sol2x;
			}
		}

		return soly;

	}

	private void intersectionImpair(Atom trait, Atom pred, Atom suc, int type) {

		/* type=0 for left and 1 for right */

		double xo1 = pred.getFx();
		double yo1 = pred.getFy();
		double xo2 = suc.getFx();
		double yo2 = suc.getFy();

		double r = 1.6;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo2;
		double b2 = -2.0 * yo2;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double a = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r * r;
		double C2 = xo2 * xo2 + yo2 * yo2 - r * r;
		double C12 = C1 - C2;
		double b, c;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);
		if (bpga) {
			b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = b * b - 4.0 * a * c;

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
		solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			soly2 = solx2;
			solx1 = -(b12 * solx1 + C12) / a12;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		double solx;
		double soly;

		if (type == 0) {
			if (solx1 <= solx2) {
				solx = solx1;
				soly = soly1;
			} else {
				solx = solx2;
				soly = soly2;
			}
		} else {
			if (solx1 >= solx2) {
				solx = solx1;
				soly = soly1;
			} else {
				solx = solx2;
				soly = soly2;
			}
		}

		trait.setFx(solx);
		trait.setFy(soly);

	}

	private void intersectionImpair2(Atom trait, Atom pred, Atom suc, int type) {

		double xo1 = pred.getFx();
		double yo1 = pred.getFy();
		double xo2 = suc.getFx();
		double yo2 = suc.getFy();

		double r = 1.6;
		double r2 = 2.77;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo2;
		double b2 = -2.0 * yo2;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double a = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r * r;
		double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
		double C12 = C1 - C2;
		double b, c;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);
		if (bpga) {
			b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = b * b - 4.0 * a * c;

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
		solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			soly2 = solx2;
			solx1 = -(b12 * solx1 + C12) / a12;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		double solx;
		double soly;

		double angle1 = angleNew(suc.getFx(), suc.getFy(), pred, solx1, soly1);
		double angle2 = angleNew(suc.getFx(), suc.getFy(), pred, solx2, soly2);

		if (type == -1) {
			if (angle1 < 0) {
				solx = solx1;
				soly = soly1;
			} else {
				solx = solx2;
				soly = soly2;
			}
		} else {
			if (angle1 > 0) {
				solx = solx1;
				soly = soly1;
			} else {
				solx = solx2;
				soly = soly2;
			}

		}

		trait.setFx(solx);
		trait.setFy(soly);

	}

	private void intersectionImpair3(Atom trait, Atom pred, Atom suc, Atom suc2) {

		double xo1 = pred.getFx();
		double yo1 = pred.getFy();
		double xo2 = suc.getFx();
		double yo2 = suc.getFy();

		double r = 1.2;
		double r2 = 2.43;

		double a1 = -2.0 * xo1;
		double b1 = -2.0 * yo1;
		double a2 = -2.0 * xo2;
		double b2 = -2.0 * yo2;
		double a12 = a1 - a2;
		double b12 = b1 - b2;
		double a = a12 * a12 + b12 * b12;

		double C1 = xo1 * xo1 + yo1 * yo1 - r * r;
		double C2 = xo2 * xo2 + yo2 * yo2 - r2 * r2;
		double C12 = C1 - C2;
		double b, c;
		boolean bpga = Math.abs(b12) >= Math.abs(a12);
		if (bpga) {
			b = 2.0 * a12 * C12 + b12 * (b1 * a2 - a1 * b2);
			c = C12 * C12 + b12 * (b1 * C2 - b2 * C1);
		} else {
			b = 2.0 * b12 * C12 + a12 * (a1 * b2 - b1 * a2);
			c = C12 * C12 + a12 * (a1 * C2 - a2 * C1);
		}
		double delta = b * b - 4.0 * a * c;

		double solx1 = 0;
		double soly1 = 0;
		double solx2 = 0;
		double soly2 = 0;

		solx1 = -0.5 * (b + 1.0 * Math.sqrt(delta)) / a;
		solx2 = -0.5 * (b + -1.0 * Math.sqrt(delta)) / a;

		if (bpga) {
			soly1 = -(a12 * solx1 + C12) / b12;
			soly2 = -(a12 * solx2 + C12) / b12;
		} else {
			soly1 = solx1;
			soly2 = solx2;
			solx1 = -(b12 * solx1 + C12) / a12;
			solx2 = -(b12 * solx2 + C12) / a12;
		}

		double solx;
		double soly;

		double angle1 = angle(suc2, suc, pred);
		double angle2 = angleNew(suc.getFx(), suc.getFy(), pred, solx1, soly1);

		if (((angle1 > 0) && (angle2 < 0)) || ((angle1 < 0) && (angle2 > 0))) {
			solx = solx1;
			soly = soly1;
		} else {
			solx = solx2;
			soly = soly2;
		}

		if ((!trait.getFlag()) || (trait.getType() == 6)) {
			trait.setFx(solx);
			trait.setFy(soly);

		}

	}

	public final double angle2(Atom a1, Atom a2) {

		double c1[] = new double[2];
		double c2[] = new double[2];
		double c3[] = new double[3];

		c1[0] = a1.getFx();
		c1[1] = a1.getFy();

		c2[0] = a2.getFx();
		c2[1] = a2.getFy();

		c3[0] = a1.getFx();
		c3[1] = 0;

		int i;
		double[] v1 = new double[2];
		double[] v2 = new double[2];
		double dot, mod1, mod2, max;
		double angle;

		for (i = 0; i < 2; i++) {
			v1[i] = c2[i] - c1[i];
			v2[i] = c3[i] - c2[i];
		}

		mod1 = v1[0] * v1[0] + v1[1] * v1[1];
		mod2 = v2[0] * v2[0] + v2[1] * v2[1];

		if (mod1 > 0.0)
			mod1 = Math.sqrt(mod1);
		else
			return (0.0F);
		if (mod2 > 0.0)
			mod2 = Math.sqrt(mod2);
		else
			return (0.0F);

		dot = v1[0] * v2[0] + v1[1] * v2[1];

		if ((max = mod1 * mod2) < dot)
			return (0.0F);

		angle = (Math.acos(dot / max));

		return angle;

	}

	public void updateAromBond() {
		for (int i = 0; i < rings.size(); i++) {
			Ring ring = (Ring) rings.get(i);
			if (ring.isAromatic()) {
				LinkedList atomId = ring.getAtoms();
				for (int j = 0; j < bonds.size(); j++) {
					Bond b = (Bond) bonds.get(j);
					Atom at1 = b.getFirst();
					if (estDedans(at1.getId(), atomId)) {
						Atom at2 = b.getSecond();
						if (estDedans(at2.getId(), atomId)) {
							b.setOrder(1);
						}
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

		if (nbre != 0) {
			res = res / nbre;
		} else {
			res = 100;
		}

		return res;
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

		for (int i = 0; i < Hbonds.size(); i++) {
			Hbond hb = (Hbond) Hbonds.get(i);
			Atom a2 = hb.getAtomRes();

			double dd = Math.sqrt((x - a2.getFx()) * (x - a2.getFx())
					+ (y - a2.getFy()) * (y - a2.getFy()));

			if (dd <= r) {

				res = res + dd;
				nbre++;
			}

		}

		for (int i = 0; i < contactList.size(); i++) {
			Contact c = (Contact) contactList.get(i);
			Atom a2 = c.getTrait();
			double dd = Math.sqrt((x - a2.getFx()) * (x - a2.getFx())
					+ (y - a2.getFy()) * (y - a2.getFy()));

			if (dd <= r) {

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

}
