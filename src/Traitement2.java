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
 This class identify the chains (liagnds, protein...) from the file2 class.
 */

import java.io.*;
import java.util.*;

public class Traitement2 {

	Molecule mol;

	public Traitement2(Molecule mol) {
		this.mol = mol;

		trait0();
		/*
		 * Assign chain class depending on heteroatom, secondary structure and
		 * residue type in chain
		 */
		trait1();
		updateConect();
		updateConectRes();
		//mol.print3();
	}

	public void trait0() {
		LinkedList l = mol.getChains();
		for (int i = 0; i < l.size(); i++) {
			Chain c = (Chain) l.get(i);
			if (c.getResidues().size() == 0)
				mol.deleteChain(c.getId());
		}
	}

	public boolean containStruc(Chain c) {
		LinkedList l = c.getResidues();
		LinkedList s = mol.getStruct();

		for (int i = 0; i < l.size(); i++) {
			Residue res = (Residue) l.get(i);
			for (int j = 0; j < s.size(); j++) {
				Structure struc = (Structure) s.get(j);
				if (struc.containRes(res))
					return true;
			}
		}
		return false;
	}

	public void updateConect() {

		for (int v = 0; v < mol.chains.size(); v++) {
			Chain c = (Chain) mol.chains.get(v);
			if (c.getType().equals("ligand")) {
				LinkedList residues = c.getResidues();
				for (int i = 0; i < residues.size(); i++) {
					Residue res = (Residue) residues.get(i);
					LinkedList atoms = res.getAtoms();
					for (int j = 0; j < atoms.size(); j++) {
						Atom at = (Atom) atoms.get(j);

						if ((res.getType().equals("aa"))
								|| (res.getType().equals("dna"))) {

							for (int t = 0; t < residues.size(); t++) {
								Residue res2 = (Residue) residues.get(t);
								LinkedList atoms2 = res2.getAtoms();
								for (int r = 0; r < atoms2.size(); r++) {

									Atom at2 = (Atom) atoms2.get(r);
									if (!at.isLinked(at2)) {
										if (at2.getId() != at.getId()) {
											double dist = at.distance(at2);
											if ((dist >= 0.4) && (dist <= 1.9)) {
												if (res.getId() == res2.getId()) {
													LinkedList l = new LinkedList();
													l.add(at2);
													mol.addAtomNeig(at, l);
												} else {
													LinkedList l = new LinkedList();
													l.add(at2);
													mol.addAtomNeig(at, l);
													LinkedList l2 = new LinkedList();
													l2.add(at);
													mol.addAtomNeig(at2, l2);
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

	public void updateConectRes() {

		for (int i = 0; i < mol.chains.size(); i++) {
			Chain chain = (Chain) mol.chains.get(i);
			if (chain.getType().equals("ligand")) {
				LinkedList residues = chain.getResidues();
				for (int j = 0; j < residues.size(); j++) {
					Residue res = (Residue) residues.get(j);
					for (int u = 0; u < mol.chains.size(); u++) {
						Chain c = (Chain) mol.chains.get(u);
						LinkedList residues2 = c.getResidues();
						for (int v = 0; v < residues2.size(); v++) {
							Residue res2 = (Residue) residues2.get(v);
							if ((chain.getName() == c.getName())
									&& (res.getId() == res2.getId())) {
							} else {
								if (res.isLinked(res2)) {

									mol.addResNeig(res.getId(),
											chain.getName(), res2.getId(), c
													.getName());
								}
							}

						}

					}
				}
			}
		}

	}

	public void trait1() {
		LinkedList l = mol.getChains();
		boolean b = false;
		for (int i = 0; i < l.size(); i++) {
			Chain c = (Chain) l.get(i);
			if (c.getName() != ' ')
				b = true;
			if (c.getId() == 0) {
			} else {
				/* Chain contains only heteroatoms */
				if (c.getName() == '?') {
					c.setType("ion");
				} else {
					if (c.getNatoms() == c.getHetero())
						c.setType("ligand2");
					else {
						if ((c.getName() == ' ') && (b == true)
								&& (!containStruc(c)))
							c.setType("ligand2");
						else {
							if (c.getNres() == c.getDNA())
								c.setType("dna");
							else {
								if (c.getNres() == c.getAA())
									c.setType("protein");
								else {
									c.setType("ligand2");
								}
							}

						}

					}
				}
			}
		}

	}

	public void trait2() {
		LinkedList l = mol.getChains();
		for (int i = 0; i < l.size(); i++) {
			Chain c = (Chain) l.get(i);
			if (c.getType().equals("protein")) {
				if (c.getNres() < 11)
					c.setType("ligand");
			}
			if (c.getType().equals("dna")) {
				if (c.getNres() < 3)
					c.setType("ligand");
			}
		}

	}

	public void trait3() {
		if (mol.getChains().size() > 1) {
			concat(mol.getChains());
		}
	}

	public void concat(LinkedList l) {
		for (int i = 0; i < l.size() - 1; i++) {
			Chain c1 = (Chain) l.get(i);

			LinkedList res1 = c1.getResidues();
			Residue r1 = (Residue) res1.getFirst();
			Residue r2 = (Residue) res1.getLast();
			for (int j = i + 1; j < l.size(); j++) {
				Chain c2 = (Chain) l.get(j);

				LinkedList res2 = c2.getResidues();
				Residue r3 = (Residue) res2.getFirst();
				Residue r4 = (Residue) res2.getLast();
				if (c1.getName() != '?') {
					if (c1.getName() == c2.getName()) {
						if (c1.getType().equals(c2.getType())) {
							if (r2.getIdInt() == r3.getIdInt() - 1) {
								mol.concat(c1.getId(), c2.getId());
								concat(mol.getChains());
								break;
							}
							if (r4.getIdInt() == r1.getIdInt() - 1) {
								mol.concat(c2.getId(), c1.getId());
								concat(mol.getChains());
								break;
							}
						}
					}
				}
			}
		}
	}

	/* check if all residues in ligands are all connected */
	public void trait4() {
		LinkedList chains = mol.getChains();
		LinkedList copy = new LinkedList();

		for (int i = 0; i < chains.size(); i++) {
			Chain chain = (Chain) chains.get(i);
			if (chain.getType().equals("ligand")) {
				copy.add(chain);
			}
		}

		for (int i = 0; i < copy.size(); i++) {
			Chain chain = (Chain) copy.get(i);

			LinkedList result = divide(chain);

			if (result.size() == 1) {
			} else {

				for (int w = 0; w < result.size(); w++) {
					LinkedList l = (LinkedList) result.get(w);
					Chain c = new Chain();
					c.setType("ligand");
					c.setId(mol.getChains().size() + 1);
					c.setName(chain.getName());
					c.setResidues(l);
					for (int u = 0; u < l.size(); u++) {
						Residue r = (Residue) l.get(u);
						r.setParent(c);
					}
					mol.addChain(c);
				}
				mol.deleteChain(chain.getId());

			}

		}
	}

	public void trait5() {

		LinkedList chains = mol.getChains();
		for (int i = 0; i < chains.size(); i++) {
			Chain chain = (Chain) chains.get(i);
			if (chain.getType().equals("ligand")) {
				LinkedList residues = chain.getResidues();
				if (residues.size() != 1) {
				} else {

					Residue res = (Residue) residues.get(0);
					if (res.type().equals("ion")) {

						res.setType();
						chain.setType("ion");
						chain.setName('?');
					} else {

						char cId = chain.getName();
						String resId = res.getId();
						LinkedList chains2 = mol.getChains();

						for (int u = 0; u < chains2.size(); u++) {

							Chain cc = (Chain) chains2.get(u);
							if (cc.getId() != chain.getId()) {
								if (cc.getName() == cId) {
									LinkedList residues2 = cc.getResidues();
									for (int v = 0; v < residues2.size(); v++) {

										Residue rr = (Residue) residues2.get(v);
										if (rr.getId().equals(resId)) {
											if (res.isLinked(rr)) {

												cc.addResidue(res);
												mol.deleteChain(chain.getId());
												break;
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

	public LinkedList divide(Chain c) {
		LinkedList result = new LinkedList();
		boolean added;
		LinkedList residues = c.getResidues();
		LinkedList l = new LinkedList();
		Residue r = (Residue) residues.get(0);
		l.add(r);
		result.add(l);
		for (int i = 1; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);
			added = false;
			for (int j = 0; j < result.size(); j++) {
				LinkedList list = (LinkedList) result.get(j);
				if (c.allConect(res, list)) {
					list.add(res);
					added = true;
					break;
				}
			}
			if (added == false) {
				LinkedList temp = new LinkedList();
				temp.add(res);
				result.add(temp);
			}
		}

		return result;
	}

}
