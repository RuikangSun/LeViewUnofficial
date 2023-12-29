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
 This class generates a molecule from a PDB file.
 */

import java.io.*;
import java.util.*;
import java.net.URL;

public class File {

	String file;
	Properties prop;

	int c;

	public File() {
	};

	public File(String file) {
		this.file = file;

	}

	public Molecule getMolecule() throws IOException {

		c = 1;
		prop = new Properties();
		try {

			URL monUrl = getClass().getResource("properties/aminoAcid.txt");

			InputStreamReader ins = new InputStreamReader(monUrl.openStream());

			BufferedReader buff = new BufferedReader(ins);

			try {
				LinkedList l = new LinkedList();
				String line;
				while ((line = buff.readLine()) != null) {
					l.add(line);
				}
				prop.setAA(l);
			} finally {
				buff.close();
			}
		} catch (IOException ioe) {
			System.out.println("Error --" + ioe.toString());
		}

		try {

			URL monUrl = getClass().getResource("properties/nucleicAcid.txt");

			InputStreamReader ins = new InputStreamReader(monUrl.openStream());

			BufferedReader buff = new BufferedReader(ins);

			try {
				LinkedList l = new LinkedList();
				String line;
				while ((line = buff.readLine()) != null) {
					l.add(line);
				}
				prop.setNA(l);
			} finally {
				buff.close();
			}
		} catch (IOException ioe) {
			System.out.println("Error --" + ioe.toString());
		}

		try {

			URL monUrl = getClass().getResource("properties/radii.txt");

			InputStreamReader ins = new InputStreamReader(monUrl.openStream());

			BufferedReader buff = new BufferedReader(ins);

			try {
				LinkedList l = new LinkedList();
				LinkedList l2 = new LinkedList();
				String line;
				while ((line = buff.readLine()) != null) {
					StringTokenizer tok = new StringTokenizer(line, ";");
					String el = tok.nextToken();
					String ra = tok.nextToken();
					double radi = Double.parseDouble(ra);
					l.add(el);
					l2.add(radi);
				}
				prop.setElements(l);
				prop.setRadii(l2);
			} finally {
				buff.close();
			}
		} catch (IOException ioe) {
			System.out.println("Error --" + ioe.toString());
		}

		try {

			URL monUrl = getClass().getResource("properties/ion.txt");

			InputStreamReader ins = new InputStreamReader(monUrl.openStream());

			BufferedReader buff = new BufferedReader(ins);

			try {
				LinkedList l = new LinkedList();
				String line;
				while ((line = buff.readLine()) != null) {
					l.add(line);
				}
				prop.setIon(l);
			} finally {
				buff.close();
			}
		} catch (IOException ioe) {
			System.out.println("Error --" + ioe.toString());
		}

		if (file.equals("")) {
			return null;
		} else {
			try {
				//String tmpo = "http://www.ebi.ac.uk/pdbe-srv/view/files/"+ file + ".pdb";
				//URL monUrl = new URL(tmpo);

				//InputStreamReader ins = new InputStreamReader(monUrl.openStream());

				//System.out.println("test");
				//BufferedReader buff = new BufferedReader(ins);
				BufferedReader buff = new BufferedReader(new FileReader(file));
				

				try {
					String line;
					Residue currentRes = null;
					Chain currentChain = null;
					Chain currentChain2 = null;
					char currentId = '*';
					int nChain = 0;
					int curId = 0;
					boolean seenENDMDL = false;
					int sss = 0;

					Molecule mol = new Molecule();
					mol.setProp(prop);
					mol.setName(file);

					while ((line = buff.readLine()) != null) {
						if (!line.equals("")) {

							char c0 = line.charAt(0);
							char c1 = line.charAt(1);
							char c2 = line.charAt(2);
							char c3 = line.charAt(3);
							char c4 = line.charAt(4);

							if (c0 == 'E' && c1 == 'N' && c2 == 'D'
									&& c3 == 'M') {
								seenENDMDL = true;
							}

							/* helix */
							if (c0 == 'H' && c1 == 'E' && c2 == 'L'
									&& c3 == 'I' && c4 == 'X') {
								char h1 = line.charAt(7);
								char h2 = line.charAt(8);
								char h3 = line.charAt(9);
								String th = "";
								if (h1 != ' ')
									th = th + h1;
								if (h2 != ' ')
									th = th + h2;
								if (h3 != ' ')
									th = th + h3;
								int helixId = Integer.parseInt(th);

								char cbeg = line.charAt(19);
								char cend = line.charAt(31);

								char hh1 = line.charAt(21);
								char hh2 = line.charAt(22);
								char hh3 = line.charAt(23);
								char hh4 = line.charAt(24);
								String tres = "";
								if (hh1 != ' ')
									tres = tres + hh1;
								if (hh2 != ' ')
									tres = tres + hh2;
								if (hh3 != ' ')
									tres = tres + hh3;
								if (hh4 != ' ')
									tres = tres + hh4;
								int begId = Integer.parseInt(tres);

								char hhh1 = line.charAt(33);
								char hhh2 = line.charAt(34);
								char hhh3 = line.charAt(35);
								char hhh4 = line.charAt(36);
								String tresh = "";
								if (hhh1 != ' ')
									tresh = tresh + hhh1;
								if (hhh2 != ' ')
									tresh = tresh + hhh2;
								if (hhh3 != ' ')
									tresh = tresh + hhh3;
								if (hhh4 != ' ')
									tresh = tresh + hhh4;
								int endId = Integer.parseInt(tresh);

								Structure struc = new Structure();
								struc.setType('H');
								struc.setId(helixId);
								struc.setBeg(begId);
								struc.setEnd(endId);
								struc.setCbeg(cbeg);
								struc.setCend(cend);
								mol.addStruct(struc);
							}

							/* sheet */
							if (c0 == 'S' && c1 == 'H' && c2 == 'E'
									&& c3 == 'E' && c4 == 'T') {
								char h1 = line.charAt(7);
								char h2 = line.charAt(8);
								char h3 = line.charAt(9);
								String th = "";
								if (h1 != ' ')
									th = th + h1;
								if (h2 != ' ')
									th = th + h2;
								if (h3 != ' ')
									th = th + h3;
								int helixId = Integer.parseInt(th);

								char cbeg = line.charAt(21);
								char cend = line.charAt(32);

								char hh1 = line.charAt(22);
								char hh2 = line.charAt(23);
								char hh3 = line.charAt(24);
								char hh4 = line.charAt(25);
								String tres = "";
								if (hh1 != ' ')
									tres = tres + hh1;
								if (hh2 != ' ')
									tres = tres + hh2;
								if (hh3 != ' ')
									tres = tres + hh3;
								if (hh4 != ' ')
									tres = tres + hh4;
								int begId = Integer.parseInt(tres);

								char hhh1 = line.charAt(33);
								char hhh2 = line.charAt(34);
								char hhh3 = line.charAt(35);
								char hhh4 = line.charAt(36);
								String tresh = "";
								if (hhh1 != ' ')
									tresh = tresh + hhh1;
								if (hhh2 != ' ')
									tresh = tresh + hhh2;
								if (hhh3 != ' ')
									tresh = tresh + hhh3;
								if (hhh4 != ' ')
									tresh = tresh + hhh4;
								int endId = Integer.parseInt(tresh);

								Structure struc = new Structure();
								struc.setType('E');
								struc.setId(helixId);
								struc.setBeg(begId);
								struc.setEnd(endId);
								struc.setCbeg(cbeg);
								struc.setCend(cend);
								mol.addStruct(struc);
							}

							/* conect record */
							if (c0 == 'C' && c1 == 'O' && c2 == 'N'
									&& c3 == 'E') {
								int debut = 0;

								LinkedList neig = new LinkedList();
								LinkedList Hneig = new LinkedList();
								if (line.charAt(10) != ' ') {
									char z1 = line.charAt(6);
									char z2 = line.charAt(7);
									char z3 = line.charAt(8);
									char z4 = line.charAt(9);
									char z5 = line.charAt(10);
									String tmp = "";
									if (z1 != ' ')
										tmp = tmp + z1;
									if (z2 != ' ')
										tmp = tmp + z2;
									if (z3 != ' ')
										tmp = tmp + z3;
									if (z4 != ' ')
										tmp = tmp + z4;
									if (z5 != ' ')
										tmp = tmp + z5;
									int zId = Integer.parseInt(tmp);
									debut = zId;
								}
								if (mol.getAtom(debut) != null) {

									Atom debb = mol.getAtom(debut);

									if (line.charAt(15) != ' ') {
										char z1 = line.charAt(11);
										char z2 = line.charAt(12);
										char z3 = line.charAt(13);
										char z4 = line.charAt(14);
										char z5 = line.charAt(15);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);

										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}
									if (line.charAt(20) != ' ') {
										char z1 = line.charAt(16);
										char z2 = line.charAt(17);
										char z3 = line.charAt(18);
										char z4 = line.charAt(19);
										char z5 = line.charAt(20);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);
										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}
									if (line.charAt(25) != ' ') {
										char z1 = line.charAt(21);
										char z2 = line.charAt(22);
										char z3 = line.charAt(23);
										char z4 = line.charAt(24);
										char z5 = line.charAt(25);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);
										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}
									if (line.charAt(30) != ' ') {
										char z1 = line.charAt(26);
										char z2 = line.charAt(27);
										char z3 = line.charAt(28);
										char z4 = line.charAt(29);
										char z5 = line.charAt(30);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);
										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}
									/* Explicit no covalent bonds */
									if (line.charAt(35) != ' ') {
										char z1 = line.charAt(31);
										char z2 = line.charAt(32);
										char z3 = line.charAt(33);
										char z4 = line.charAt(34);
										char z5 = line.charAt(35);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);
										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}
									if (line.charAt(40) != ' ') {
										char z1 = line.charAt(36);
										char z2 = line.charAt(37);
										char z3 = line.charAt(38);
										char z4 = line.charAt(39);
										char z5 = line.charAt(40);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);
										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}
									if (line.charAt(45) != ' ') {
										char z1 = line.charAt(41);
										char z2 = line.charAt(42);
										char z3 = line.charAt(43);
										char z4 = line.charAt(44);
										char z5 = line.charAt(45);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);
										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}
									if (line.charAt(50) != ' ') {
										char z1 = line.charAt(46);
										char z2 = line.charAt(47);
										char z3 = line.charAt(48);
										char z4 = line.charAt(49);
										char z5 = line.charAt(50);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);
										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}
									if (line.charAt(55) != ' ') {
										char z1 = line.charAt(51);
										char z2 = line.charAt(52);
										char z3 = line.charAt(53);
										char z4 = line.charAt(54);
										char z5 = line.charAt(55);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);
										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}
									if (line.charAt(60) != ' ') {
										char z1 = line.charAt(56);
										char z2 = line.charAt(57);
										char z3 = line.charAt(58);
										char z4 = line.charAt(59);
										char z5 = line.charAt(60);
										String tmp = "";
										if (z1 != ' ')
											tmp = tmp + z1;
										if (z2 != ' ')
											tmp = tmp + z2;
										if (z3 != ' ')
											tmp = tmp + z3;
										if (z4 != ' ')
											tmp = tmp + z4;
										if (z5 != ' ')
											tmp = tmp + z5;
										int zId = Integer.parseInt(tmp);
										Atom atom = mol.getAtom(zId);
										if (atom != null) {
											if (atom.getParent().getId()
													.equals(
															debb.getParent()
																	.getId())) {
												if (atom != null) {
													neig.add(atom);
												}
											} else {
												if ((atom.getParent()
														.getIdInt() - 1 == debb
														.getParent().getIdInt())
														|| (atom.getParent()
																.getIdInt() + 1 == debb
																.getParent()
																.getIdInt())) {
													if (atom != null) {
														neig.add(atom);
													}
												} else {
													if (atom != null)
														Hneig.add(atom);
												}
											}
										}
									}

									if (Hneig.size() != 0) {
										mol.addAtomHneig(mol.getAtom(debut),
												Hneig);
									}

									mol.addAtomNeig(mol.getAtom(debut), neig);
								}

							}

							if (c0 == 'T' && c1 == 'E' && c2 == 'R') {
								currentId = '*';
							}

							/* ATOM OR HETATM */
							if ((seenENDMDL == false)
									&& ((c0 == 'A' && c1 == 'T' && c2 == 'O' && c3 == 'M') || (c0 == 'H'
											&& c1 == 'E' && c2 == 'T' && c3 == 'A'))) {

								char alt = line.charAt(16);
								if ((alt == ' ') || (alt == 'A')) {

									char ca = line.charAt(17);
									char cb = line.charAt(18);
									char cc = line.charAt(19);
									String residueName = "";
									if (ca != ' ')
										residueName = residueName + ca;
									if (cb != ' ')
										residueName = residueName + cb;
									if (cc != ' ')
										residueName = residueName + cc;

									char caa = line.charAt(76);
									char cab = line.charAt(77);
									String atomElement = "";
									if (caa != ' ')
										atomElement = atomElement + caa;
									if (cab != ' ')
										atomElement = atomElement + cab;

									/* ignore water and hydrogens */
									if ((!isWater(residueName))
											&& (!atomElement.equals("H"))) {

										/* get Residue Id */
										String temp = "";
										if (line.charAt(22) != ' ')
											temp = temp + line.charAt(22);
										if (line.charAt(23) != ' ')
											temp = temp + line.charAt(23);
										if (line.charAt(24) != ' ')
											temp = temp + line.charAt(24);
										if (line.charAt(25) != ' ')
											temp = temp + line.charAt(25);
										if (line.charAt(26) != ' ')
											temp = temp + line.charAt(26);
										String residueId = temp;

										char insertionCode = line.charAt(26);
										char chainId = line.charAt(21);

										if (chainId == currentId) {
										} else {
											nChain = nChain + 1;
											Chain chain = new Chain();
											chain.setName(chainId);
											chain.setId(nChain);
											chain.setParent(mol);
											mol.addChain(chain);
											currentId = chainId;
											currentChain = chain;
											curId = nChain;

										}

										if (currentChain.containRes(residueId)) {
										} else {
											Residue res = new Residue();
											res.setName(residueName);
											res.setId(residueId);
											if (res.isIon(currentChain)) {
												nChain = nChain + 1;
												Chain chain = new Chain();
												chain.setParent(mol);
												chain.setName('?');
												chain.setId(nChain);
												mol.addChain(chain);

												currentChain2 = chain;
												res.setParent(chain);
												mol.addResidue(nChain, res);

												currentRes = res;

											} else {
												res.setParent(currentChain);
												mol.addResidue(curId, res);
												currentRes = res;
											}
										}

										// atom name
										char c12 = line.charAt(12);
										char c13 = line.charAt(13);
										char c14 = line.charAt(14);
										char c15 = line.charAt(15);
										String atomLabel = "";
										if (c12 != ' ')
											atomLabel = atomLabel + c12;
										if (c13 != ' ')
											atomLabel = atomLabel + c13;
										if (c14 != ' ')
											atomLabel = atomLabel + c14;
										if (c15 != ' ')
											atomLabel = atomLabel + c15;

										// atom id
										char cba = line.charAt(6);
										char cbb = line.charAt(7);
										char cbc = line.charAt(8);
										char cbd = line.charAt(9);
										char cbe = line.charAt(10);
										String tempId = "";
										if (cba != ' ')
											tempId = tempId + cba;
										if (cbb != ' ')
											tempId = tempId + cbb;
										if (cbc != ' ')
											tempId = tempId + cbc;
										if (cbd != ' ')
											tempId = tempId + cbd;
										if (cbe != ' ')
											tempId = tempId + cbe;
										int atomId = Integer.parseInt(tempId);

										// coordinates
										char x1 = line.charAt(30);
										char x2 = line.charAt(31);
										char x3 = line.charAt(32);
										char x4 = line.charAt(33);
										char x5 = line.charAt(34);
										char x6 = line.charAt(35);
										char x7 = line.charAt(36);
										char x8 = line.charAt(37);
										String tempX = "";
										if (x1 != ' ')
											tempX = tempX + x1;
										if (x2 != ' ')
											tempX = tempX + x2;
										if (x3 != ' ')
											tempX = tempX + x3;
										if (x4 != ' ')
											tempX = tempX + x4;
										if (x5 != ' ')
											tempX = tempX + x5;
										if (x6 != ' ')
											tempX = tempX + x6;
										if (x7 != ' ')
											tempX = tempX + x7;
										if (x8 != ' ')
											tempX = tempX + x8;
										double x = Double.parseDouble(tempX);

										char y1 = line.charAt(38);
										char y2 = line.charAt(39);
										char y3 = line.charAt(40);
										char y4 = line.charAt(41);
										char y5 = line.charAt(42);
										char y6 = line.charAt(43);
										char y7 = line.charAt(44);
										char y8 = line.charAt(45);
										String tempY = "";
										if (y1 != ' ')
											tempY = tempY + y1;
										if (y2 != ' ')
											tempY = tempY + y2;
										if (y3 != ' ')
											tempY = tempY + y3;
										if (y4 != ' ')
											tempY = tempY + y4;
										if (y5 != ' ')
											tempY = tempY + y5;
										if (y6 != ' ')
											tempY = tempY + y6;
										if (y7 != ' ')
											tempY = tempY + y7;
										if (y8 != ' ')
											tempY = tempY + y8;
										double y = Double.parseDouble(tempY);

										char z1 = line.charAt(46);
										char z2 = line.charAt(47);
										char z3 = line.charAt(48);
										char z4 = line.charAt(49);
										char z5 = line.charAt(50);
										char z6 = line.charAt(51);
										char z7 = line.charAt(52);
										char z8 = line.charAt(53);
										String tempZ = "";
										if (z1 != ' ')
											tempZ = tempZ + z1;
										if (z2 != ' ')
											tempZ = tempZ + z2;
										if (z3 != ' ')
											tempZ = tempZ + z3;
										if (z4 != ' ')
											tempZ = tempZ + z4;
										if (z5 != ' ')
											tempZ = tempZ + z5;
										if (z6 != ' ')
											tempZ = tempZ + z6;
										if (z7 != ' ')
											tempZ = tempZ + z7;
										if (z8 != ' ')
											tempZ = tempZ + z8;
										double z = Double.parseDouble(tempZ);

										Atom atom = new Atom();
										atom.setName(atomLabel);
										atom.setId(atomId);
										atom.setElement(atomElement);
										atom.setCoord(x, y, z);
										atom.setParent(currentRes);
										if (c0 == 'A')
											atom.setHetero(false);
										else
											atom.setHetero(true);
										if (currentRes.isIon(currentChain)) {
											mol.addAtom(currentChain2,
													currentRes, atom);
										} else {
											mol.addAtom(currentChain,
													currentRes, atom);
										}
										currentRes.setType();

									} else {

										if (isWater(residueName)) {
											/* get Residue Id */
											String temp = "";
											if (line.charAt(22) != ' ')
												temp = temp + line.charAt(22);
											if (line.charAt(23) != ' ')
												temp = temp + line.charAt(23);
											if (line.charAt(24) != ' ')
												temp = temp + line.charAt(24);
											if (line.charAt(25) != ' ')
												temp = temp + line.charAt(25);
											if (line.charAt(26) != ' ')
												temp = temp + line.charAt(26);
											String residueId = temp;

											char insertionCode = line
													.charAt(26);
											char chainId = line.charAt(21);

											// atom name
											char c12 = line.charAt(12);
											char c13 = line.charAt(13);
											char c14 = line.charAt(14);
											char c15 = line.charAt(15);
											String atomLabel = "";
											if (c12 != ' ')
												atomLabel = atomLabel + c12;
											if (c13 != ' ')
												atomLabel = atomLabel + c13;
											if (c14 != ' ')
												atomLabel = atomLabel + c14;
											if (c15 != ' ')
												atomLabel = atomLabel + c15;

											// atom id
											char cba = line.charAt(6);
											char cbb = line.charAt(7);
											char cbc = line.charAt(8);
											char cbd = line.charAt(9);
											char cbe = line.charAt(10);
											String tempId = "";
											if (cba != ' ')
												tempId = tempId + cba;
											if (cbb != ' ')
												tempId = tempId + cbb;
											if (cbc != ' ')
												tempId = tempId + cbc;
											if (cbd != ' ')
												tempId = tempId + cbd;
											if (cbe != ' ')
												tempId = tempId + cbe;
											int atomId = Integer
													.parseInt(tempId);

											// coordinates
											char x1 = line.charAt(30);
											char x2 = line.charAt(31);
											char x3 = line.charAt(32);
											char x4 = line.charAt(33);
											char x5 = line.charAt(34);
											char x6 = line.charAt(35);
											char x7 = line.charAt(36);
											char x8 = line.charAt(37);
											String tempX = "";
											if (x1 != ' ')
												tempX = tempX + x1;
											if (x2 != ' ')
												tempX = tempX + x2;
											if (x3 != ' ')
												tempX = tempX + x3;
											if (x4 != ' ')
												tempX = tempX + x4;
											if (x5 != ' ')
												tempX = tempX + x5;
											if (x6 != ' ')
												tempX = tempX + x6;
											if (x7 != ' ')
												tempX = tempX + x7;
											if (x8 != ' ')
												tempX = tempX + x8;
											double x = Double
													.parseDouble(tempX);

											char y1 = line.charAt(38);
											char y2 = line.charAt(39);
											char y3 = line.charAt(40);
											char y4 = line.charAt(41);
											char y5 = line.charAt(42);
											char y6 = line.charAt(43);
											char y7 = line.charAt(44);
											char y8 = line.charAt(45);
											String tempY = "";
											if (y1 != ' ')
												tempY = tempY + y1;
											if (y2 != ' ')
												tempY = tempY + y2;
											if (y3 != ' ')
												tempY = tempY + y3;
											if (y4 != ' ')
												tempY = tempY + y4;
											if (y5 != ' ')
												tempY = tempY + y5;
											if (y6 != ' ')
												tempY = tempY + y6;
											if (y7 != ' ')
												tempY = tempY + y7;
											if (y8 != ' ')
												tempY = tempY + y8;
											double y = Double
													.parseDouble(tempY);

											char z1 = line.charAt(46);
											char z2 = line.charAt(47);
											char z3 = line.charAt(48);
											char z4 = line.charAt(49);
											char z5 = line.charAt(50);
											char z6 = line.charAt(51);
											char z7 = line.charAt(52);
											char z8 = line.charAt(53);
											String tempZ = "";
											if (z1 != ' ')
												tempZ = tempZ + z1;
											if (z2 != ' ')
												tempZ = tempZ + z2;
											if (z3 != ' ')
												tempZ = tempZ + z3;
											if (z4 != ' ')
												tempZ = tempZ + z4;
											if (z5 != ' ')
												tempZ = tempZ + z5;
											if (z6 != ' ')
												tempZ = tempZ + z6;
											if (z7 != ' ')
												tempZ = tempZ + z7;
											if (z8 != ' ')
												tempZ = tempZ + z8;
											double z = Double
													.parseDouble(tempZ);

											Atom atom = new Atom();
											atom.setName(atomLabel);
											atom.setId(atomId);
											atom.setElement(atomElement);
											atom.setCoord(x, y, z);
											atom.setWaterId(c);

											mol.addWater(atom);
											c++;

										}

									}

								}
							}

						}

					}

					return mol;
				} finally {
					buff.close();
				}
			} catch (IOException ioe) {

				System.out.println("Error --" + ioe.toString());
				return null;
			}
		}

	}

	public boolean isWater(String s) {
		if ((s.equals("HOH")) || (s.equals("H2O")) || (s.equals("WAT")))
			return true;
		else
			return false;
	}

}
