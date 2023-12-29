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
 This class defines the panel for the ligand frame.
 */

import java.awt.Color;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.*;

import java.awt.BasicStroke;
 
import javax.swing.JPanel;
import javax.swing.*;

import java.text.*; 
import javax.swing.SwingUtilities;


import java.io.*;
import java.util.*;

class MyPanel extends JPanel implements MouseListener, MouseMotionListener {

	Ligand ligand;
	LinkedList atoms;
	double[] viewSize;
	LinkedList residues;
	LinkedList bonds;
	LinkedList explicit;
	LinkedList rings;
	Atom[] view;
	LinkedList Hbonds;
	LinkedList contactList;
	MyFrame frame;
	TopMenu panel;
	static int x, y;

	LinkedList resExplicit;

	LinkedList posName;
	LinkedList posName2;

	double cutoffHB;
	double cutoffC;
	double valueHB;
	double valueC;

	int colorMode;
	int mode;
	Contact traitContact;
	Hbond traitHbond;
	Explicit traitExplicit;
	Atom atomExplicit;
	int traitPosName;
	int traitPosName2;
	Atom traitLigand;

	int lockedMode;

	LinkedList waterMediated;
	LinkedList currentWaterMediated;
	Path traitPath;

	Color ligandColor;
	Color resExplicitColor;
	static Color hbColor;
	Color wmColor;
	Color resColor;
	Color nameColor;
	Color name2Color;

	int hbStyle;
	int labelMode;
	int atomColorMode;
	int displayDist;


	Atom selectedAtom=null;

	public MyPanel() {
	}

	public MyPanel(MyFrame f, TopMenu t) {
		frame = f;
		panel = t;
		addMouseListener(this);
		addMouseMotionListener(this);
		mode = 0;
		lockedMode = 0;
	}

	double scale;

	public JPanel affiche(Ligand ligand) {
		this.ligand = ligand;
		//ligand.print();
		waterMediated = ligand.getMediatedPath();
		atoms = ligand.getAtoms();
		bonds = ligand.getBonds();
		Hbonds = ligand.getHbond();
		resExplicit = ligand.getResExpli();
		explicit = ligand.getExplicit();
		residues = ligand.getResidues();
		contactList = ligand.getContacts();
		rings = ligand.getRing();
		traitPosName = -1;
		traitPosName2 = -1;
		posName = ligand.returnPoseName();
		posName2 = ligand.returnPoseName2();
		currentWaterMediated = new LinkedList();

		init();

		setScale();
		cutoffHB = 3.3;
		cutoffC = 4.0;
		colorMode = 1;

		ligandColor = Color.BLACK;
		resExplicitColor = Color.GRAY;
		hbColor = Color.BLUE;
		wmColor = Color.PINK;
		resColor = Color.WHITE;

		nameColor = Color.RED;
		name2Color = Color.GRAY;

		int hbStyle = 1;
		labelMode = 1;
		atomColorMode = 0;
		displayDist = 0;

		//this.setBounds(10, 10, 300, 300);
		this.setBackground(Color.WHITE);
		return this;
	}

	public boolean isInsideEx(Explicit e, LinkedList l) {
		Residue r1 = e.getResidue();
		for (int i = 0; i < l.size(); i++) {
			Explicit ex = (Explicit) l.get(i);
			Residue r2 = ex.getResidue();
			if (r1.getId() == r2.getId())
				return true;
		}
		return false;
	}

	public void init() {
		LinkedList temp = new LinkedList();
		/* copying atoms */
		for (int i = 0; i < atoms.size(); i++) {
			Atom a1 = (Atom) atoms.get(i);
			temp.add(a1.getCopy());
		}
		atoms = temp;
		temp = new LinkedList();
		/* copying bonds */
		for (int i = 0; i < bonds.size(); i++) {
			Bond b = (Bond) bonds.get(i);
			Bond b1 = new Bond(getCorres(b.getFirst()),
					getCorres(b.getSecond()), b.getOrder());
			temp.add(b1);
		}
		bonds = temp;
		/* copying Hbonds */
		LinkedList tempHB = new LinkedList();
		temp = new LinkedList();
		for (int i = 0; i < Hbonds.size(); i++) {
			Hbond b = (Hbond) Hbonds.get(i);
			Atom don = b.getDon();
			Atom acc = b.getAcc();
			if (getCorres(don) != null) {

				if (getCorres(acc, tempHB) == null) {
					Atom ttt = acc.getCopy();
					Hbond b1 = new Hbond(getCorres(don), ttt, b.getParent(), b
							.getDistance());
					tempHB.add(ttt);
					temp.add(b1);
				} else {
					Hbond b1 = new Hbond(getCorres(don),
							getCorres(acc, tempHB), b.getParent(), b
									.getDistance());
					temp.add(b1);
				}
			} else {
				if (getCorres(don, tempHB) == null) {
					Atom ttt = don.getCopy();
					Hbond b1 = new Hbond(ttt, getCorres(acc), b.getParent(), b
							.getDistance());
					tempHB.add(ttt);
					temp.add(b1);
				} else {
					Hbond b1 = new Hbond(getCorres(don, tempHB),
							getCorres(acc), b.getParent(), b.getDistance());
					temp.add(b1);
				}
			}
		}
		Hbonds = temp;
		/* copying resExplicit */
		temp = new LinkedList();
		for (int i = 0; i < resExplicit.size(); i++) {
			Explicit e = (Explicit) resExplicit.get(i);
			Explicit e2 = e.getCopy();
			e2.setAtLigand(getCorres(e.getAtLigand()));
			if (!isInsideEx(e2, temp)) {
				temp.add(e2);
			}
		}
		resExplicit = temp;
		temp = new LinkedList();

		for (int i = 0; i < explicit.size(); i++) {
			Bond b = (Bond) explicit.get(i);

			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();
			if (getCorres(a1) != null) {
				Bond b1 = new Bond(getCorres(a1), getCorres2(a2), b.getOrder());
				temp.add(b1);
			} else {
				Bond b1 = new Bond(getCorres2(a1), getCorres(a2), b.getOrder());
				temp.add(b1);
			}
		}
		explicit = temp;

		temp = new LinkedList();

		for (int i = 0; i < waterMediated.size(); i++) {
			Path p = (Path) waterMediated.get(i);
			p.setDep(getCorres(p.getDep()));
		}
		temp = new LinkedList();
		for (int i = 0; i < contactList.size(); i++) {
			Contact c = (Contact) contactList.get(i);
			Contact c2 = new Contact(c.getResidue(), c.getTrait().getCopy(), c
					.getRef().getCopy(), c.getDistance() / 2);
			temp.add(c2);
		}
		contactList = temp;

	}

	public Atom getCorres(Atom at, LinkedList l) {
		for (int i = 0; i < l.size(); i++) {
			Atom a1 = (Atom) l.get(i);
			if (at.getId() == a1.getId())
				return a1;
		}
		return null;

	}

	public void deleteHbond(Hbond h){

	for(int i=0;i<Hbonds.size();i++){
		Hbond hh=(Hbond)Hbonds.get(i);
		if(h.equals(hh)){
		Hbonds.remove(i);
		return;
		}

	}

	}


	public void deleteContact(Contact c){

	for(int i=0;i<contactList.size();i++){
		Contact cc=(Contact)contactList.get(i);
		if(c.equals(cc)){
		contactList.remove(i);
		return;
		}

	}

	}

	public Atom getCorres(Atom at) {

		for (int i = 0; i < atoms.size(); i++) {
			Atom a1 = (Atom) atoms.get(i);
			if (at.getId() == a1.getId())
				return a1;
		}
		return null;
	}

	public Atom getCorres2(Atom at) {

		for (int i = 0; i < resExplicit.size(); i++) {
			Explicit e = (Explicit) resExplicit.get(i);
			LinkedList l = e.getAtoms();
			for (int j = 0; j < l.size(); j++) {
				Atom a1 = (Atom) l.get(j);
				if (at.getId() == a1.getId())
					return a1;
			}
		}
		return null;
	}

	public void setLigandColor(Color c) {
		ligandColor = c;
	}

	public void setResExplicitColor(Color c) {
		resExplicitColor = c;
	}

	public void setHbondsColor(Color c) {
		hbColor = c;
	}

	public void setWmColor(Color c) {
		wmColor = c;
	}

	public void setResColor(Color c) {
		resColor = c;
	}

	public void setNameColor(Color c) {
		nameColor = c;
	}

	public void setName2Color(Color c) {
		name2Color = c;
	}

	public void setHbStyle(int i) {
		hbStyle = i;
	}

	public void setLabelMode(int i) {
		labelMode = i;
	}

	public int getLabelMode() {
		return labelMode;
	}

	public void setAtomColorMode(int i) {
		atomColorMode = i;
	}

	public int getAtomColorMode() {
		return atomColorMode;
	}

	public void setDisplayDist(int i) {
		displayDist = i;
	}

	public int getDisplayDist() {
		return displayDist;
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		setScale();
		drawScale3();
		paintAtoms(g);
	}

	public void paintAtoms(Graphics g) {

		Graphics2D g2d = (Graphics2D) g;
		g2d.setStroke(new BasicStroke(2.0f));
		/* bonds */
		for (int i = 0; i < bonds.size(); i++) {
			Bond b = (Bond) bonds.get(i);
			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();

			// a1.setXg(getCx(a1.getFx()));
			// a1.setYg(getCy(a1.getFy()));
			// a2.setXg(getCx(a2.getFx()));
			// a2.setYg(getCy(a2.getFy()));

			int x1 = getCx(a1.getFx());
			int y1 = getCy(a1.getFy());
			int x2 = getCx(a2.getFx());
			int y2 = getCy(a2.getFy());

			if (ligand.getSelected()) {
				g.setColor(Color.RED);
			} else {
				g.setColor(ligandColor);
			}

			if (b.getOrder() == 1) {
				g.drawLine(x1 + 10, y1 + 5, x2 + 10, y2 + 5);
			}
			if (b.getOrder() == 2) {
				LinkedList l1 = getDouble(a1, a2);
				LinkedList l2 = getDouble(a2, a1);

				int X1 = 0;
				int Y1 = 0;
				int X2 = 0;
				int Y2 = 0;
				int X3 = 0;
				int Y3 = 0;
				int X4 = 0;
				int Y4 = 0;

				double xx1 = (Double) l1.get(0);
				double xx2 = (Double) l1.get(2);
				double yy1 = (Double) l1.get(1);
				double yy2 = (Double) l1.get(3);

				if (xx1 != xx2) {
					if (xx1 < xx2) {
						X1 = getCx(xx1);
						Y1 = getCy((Double) l1.get(1));
						X2 = getCx(xx2);
						Y2 = getCy((Double) l1.get(3));
					} else {
						X1 = getCx(xx2);
						Y1 = getCy((Double) l1.get(3));
						X2 = getCx(xx1);
						Y2 = getCy((Double) l1.get(1));
					}
				} else {
					if (yy1 < yy2) {
						X1 = getCx(xx1);
						Y1 = getCy((Double) l1.get(1));
						X2 = getCx(xx2);
						Y2 = getCy((Double) l1.get(3));
					} else {
						X1 = getCx(xx2);
						Y1 = getCy((Double) l1.get(3));
						X2 = getCx(xx1);
						Y2 = getCy((Double) l1.get(1));
					}
				}

				double xxx1 = (Double) l2.get(0);
				double xxx2 = (Double) l2.get(2);
				double yyy1 = (Double) l2.get(1);
				double yyy2 = (Double) l2.get(3);

				if (xxx1 != xxx2) {
					if (xxx1 < xxx2) {
						X3 = getCx(xxx1);
						Y3 = getCy((Double) l2.get(1));
						X4 = getCx(xxx2);
						Y4 = getCy((Double) l2.get(3));
					} else {
						X3 = getCx(xxx2);
						Y3 = getCy((Double) l2.get(3));
						X4 = getCx(xxx1);
						Y4 = getCy((Double) l2.get(1));
					}
				} else {
					if (yyy1 < yyy2) {
						X3 = getCx(xxx1);
						Y3 = getCy((Double) l2.get(1));
						X4 = getCx(xxx2);
						Y4 = getCy((Double) l2.get(3));
					} else {
						X3 = getCx(xxx2);
						Y3 = getCy((Double) l2.get(3));
						X4 = getCx(xxx1);
						Y4 = getCy((Double) l2.get(1));
					}
				}

				g.drawLine(X1 + 10, Y1 + 5, X3 + 10, Y3 + 5);
				g.drawLine(X2 + 10, Y2 + 5, X4 + 10, Y4 + 5);

			}
			if (b.getOrder() == 3) {

				g.drawLine(x1 + 10, y1 + 5, x2 + 10, y2 + 5);

				LinkedList l1 = getDouble(a1, a2);
				LinkedList l2 = getDouble(a2, a1);

				int X1 = 0;
				int Y1 = 0;
				int X2 = 0;
				int Y2 = 0;
				int X3 = 0;
				int Y3 = 0;
				int X4 = 0;
				int Y4 = 0;

				double xx1 = (Double) l1.get(0);
				double xx2 = (Double) l1.get(2);
				double yy1 = (Double) l1.get(1);
				double yy2 = (Double) l1.get(3);

				if (xx1 != xx2) {
					if (xx1 < xx2) {
						X1 = getCx(xx1);
						Y1 = getCy((Double) l1.get(1));
						X2 = getCx(xx2);
						Y2 = getCy((Double) l1.get(3));
					} else {
						X1 = getCx(xx2);
						Y1 = getCy((Double) l1.get(3));
						X2 = getCx(xx1);
						Y2 = getCy((Double) l1.get(1));
					}
				} else {
					if (yy1 < yy2) {
						X1 = getCx(xx1);
						Y1 = getCy((Double) l1.get(1));
						X2 = getCx(xx2);
						Y2 = getCy((Double) l1.get(3));
					} else {
						X1 = getCx(xx2);
						Y1 = getCy((Double) l1.get(3));
						X2 = getCx(xx1);
						Y2 = getCy((Double) l1.get(1));
					}
				}

				double xxx1 = (Double) l2.get(0);
				double xxx2 = (Double) l2.get(2);
				double yyy1 = (Double) l2.get(1);
				double yyy2 = (Double) l2.get(3);

				if (xxx1 != xxx2) {
					if (xxx1 < xxx2) {
						X3 = getCx(xxx1);
						Y3 = getCy((Double) l2.get(1));
						X4 = getCx(xxx2);
						Y4 = getCy((Double) l2.get(3));
					} else {
						X3 = getCx(xxx2);
						Y3 = getCy((Double) l2.get(3));
						X4 = getCx(xxx1);
						Y4 = getCy((Double) l2.get(1));
					}
				} else {
					if (yyy1 < yyy2) {
						X3 = getCx(xxx1);
						Y3 = getCy((Double) l2.get(1));
						X4 = getCx(xxx2);
						Y4 = getCy((Double) l2.get(3));
					} else {
						X3 = getCx(xxx2);
						Y3 = getCy((Double) l2.get(3));
						X4 = getCx(xxx1);
						Y4 = getCy((Double) l2.get(1));
					}
				}

				g.drawLine(X1 + 10, Y1 + 5, X3 + 10, Y3 + 5);
				g.drawLine(X2 + 10, Y2 + 5, X4 + 10, Y4 + 5);

			}
		}

		g2d.setStroke(new BasicStroke(1.0f));

		/* explicit bonds */
		for (int i = 0; i < explicit.size(); i++) {
			Bond b = (Bond) explicit.get(i);
			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();

			// a1.setXg(getCx(a1.getFx()));
			// a1.setYg(getCy(a1.getFy()));
			// a2.setXg(getCx(a2.getFx()));
			// a2.setYg(getCy(a2.getFy()));

			int x1 = getCx(a1.getFx());
			int y1 = getCy(a1.getFy());
			int x2 = getCx(a2.getFx());
			int y2 = getCy(a2.getFy());

			g.setColor(Color.RED);

			if (b.getOrder() == 1) {
				g.drawLine(x1 + 10, y1 + 5, x2 + 10, y2 + 5);
			}
			if (b.getOrder() == 2) {
				g.drawLine(x1 + 12, y1 + 7, x2 + 12, y2 + 7);
				g.drawLine(x1 + 8, y1 + 3, x2 + 8, y2 + 3);
			}
		}

		/* atoms */
		if (labelMode == 0) {
			for (int i = 0; i < atoms.size(); i++) {

				Atom atom = (Atom) atoms.get(i);

				// atom.setXg(getCx(atom.getFx()));
				// atom.setYg(getCy(atom.getFy()));

				int x = getCx(atom.getFx());
				int y = getCy(atom.getFy());

				g.setColor(Color.white);
				g.fillRect(x, y, 20, 10);
				g.setColor(Color.white);
				g.drawRect(x, y, 20, 10);

				if (ligand.getSelected()) {
					g.setColor(Color.RED);
				} else {
					if (atomColorMode == 1) {
						g.setColor(ligandColor);
					} else {
						g.setColor(getAtomColor(atom.getElement()));
					}
				}

				int leng;
				if (atom.getName().length() == 1) {
					leng = 6;
				} else
					leng = 1;
				g.drawString(atom.getName(), x + leng, y + 10);

			}

		} else {

			for (int i = 0; i < atoms.size(); i++) {

				Atom atom = (Atom) atoms.get(i);

				if (!atom.getElement().equals("C")) {

					// atom.setXg(getCx(atom.getFx()));
					// atom.setYg(getCy(atom.getFy()));

					int x = getCx(atom.getFx());
					int y = getCy(atom.getFy());

					g.setColor(Color.white);
					g.fillRect(x + 5, y, 10, 10);
					g.setColor(Color.white);
					g.drawRect(x + 5, y, 10, 10);

					if (ligand.getSelected()) {
						g.setColor(Color.RED);
					} else {
						if (atomColorMode == 1) {
							g.setColor(ligandColor);
						} else {
							g.setColor(getAtomColor(atom.getElement()));
						}
					}

					g.drawString(atom.getElement(), x + 5, y + 10);

				}
				/* Carbon */
				else {
					// atom.setXg(getCx(atom.getFx()));
					// atom.setYg(getCy(atom.getFy()));

					int x = getCx(atom.getFx());
					int y = getCy(atom.getFy());

					if (ligand.getSelected()) {
						g.setColor(Color.RED);
					} else {
						if (atomColorMode == 1) {
							g.setColor(ligandColor);
						} else {
							g.setColor(getAtomColor(atom.getElement()));
						}
					}

					g.drawOval(x + 7, y + 2, 5, 5);
					g.fillOval(x + 7, y + 2, 5, 5);

				}

			}

		}

		if (ligand.getSelected()) {
			g.setColor(Color.RED);
		} else {
			g.setColor(ligandColor);
		}

		/* aromatic */
		g2d.setStroke(new BasicStroke(2.0f));
		for (int i = 0; i < rings.size(); i++) {
			Ring r = (Ring) rings.get(i);
			if (r.isAromatic()) {
				LinkedList att = r.getAtoms();
				int id = (Integer) att.get(0);

				Atom aaa = getAtom(id);

				double xo = r.getCenterX();
				double yo = r.getCenterY();

				int x1 = getCx(xo);
				int y1 = getCy(yo);

				int x2 = getCx(aaa.getFx());
				int y2 = getCy(aaa.getFy());

				double d = Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2)
						* (y1 - y2));
				int dd = (int) Math.rint(d);
				

				g.drawOval(x1 + 10 - (dd / 2), y1 + 5 - (dd / 2), dd, dd);

			}
		}
		g2d.setStroke(new BasicStroke(1.0f));

		for (int i = 0; i < contactList.size(); i++) {
			Contact con = (Contact) contactList.get(i);
			paintContact(g, con);
		}

		for (int i = 0; i < resExplicit.size(); i++) {
			Explicit e = (Explicit) resExplicit.get(i);
			paintEx(g, e);

		}

		for (int i = 0; i < Hbonds.size(); i++) {
			Hbond h = (Hbond) Hbonds.get(i);
			paintHbond(g, h);

		}

		paintName(g);

		paintName2(g);

		for (int i = 0; i < currentWaterMediated.size(); i++) {
			Path p = (Path) currentWaterMediated.get(i);

			paintPath(g, p);

		}

	}

	public LinkedList getDouble(Atom a1, Atom a2) {
		LinkedList res = new LinkedList();

		double x1 = a1.getFx();
		double y1 = a1.getFy();
		double x2 = a2.getFx();
		double y2 = a2.getFy();

		double aa = y2 - y1;
		double bb = x1 - x2;
		double cc = y1 * x2 - x1 * y2;

		double a = bb;
		double b = -aa;
		double tt = a * x1 + b * y1;
		double c = -tt;

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double xo = x1;
		double yo = y1;
		double r = 0.1;

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

		res.add(sol1x);
		res.add(sol1y);
		res.add(sol2x);
		res.add(sol2y);

		return res;

	}

	public Color getAtomColor(String e) {

		if (e.equals("C"))
			return Color.BLACK;
		if (e.equals("O"))
			return Color.RED;
		if (e.equals("N"))
			return Color.BLUE;
		if (e.equals("S"))
			return new Color(255, 204, 0);
		if (e.equals("Cl"))
			return Color.GREEN;

		return Color.GRAY;

	}

	public Atom getAtom(int id) {
		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			if (at.getId() == id)
				return at;
		}
		return null;
	}

	public void setCutoffHB(double cut) {
		if (lockedMode == 0) {
			ligand = new Ligand(ligand.getParent(), cutoffC, cut);
			atoms = ligand.getAtoms();
			bonds = ligand.getBonds();
			Hbonds = ligand.getHbond();
			resExplicit = ligand.getResExpli();
			explicit = ligand.getExplicit();
			residues = ligand.getResidues();
			contactList = ligand.getContacts();
			rings = ligand.getRing();
			posName = ligand.returnPoseName();
			posName2 = ligand.returnPoseName2();
			init();
			setScale();
		}
		cutoffHB = cut;
	}

	public LinkedList getWaterMediated() {
		return waterMediated;
	}

	public void setLockedMode(int i) {
		lockedMode = i;
		valueHB = cutoffHB;
		valueC = cutoffC;
	}

	public int getLockedMode() {
		return lockedMode;
	}

	public double getValueHB() {
		return valueHB;
	}

	public double getValueC() {
		return valueC;
	}

	public void setValueHB(double d) {
		valueHB = d;
	}

	public void setValueC(double d) {
		valueC = d;
	}

	public void paintName(Graphics g) {

		for (int i = 0; i < residues.size(); i++) {

			Residue r = (Residue) residues.get(i);
			Chain c = r.getParent();

			String temp = r.getName() + " " + r.getId() + "(" + c.getName()
					+ ")";
			double[] tab = (double[]) posName.get(i);

			int x = getCx(tab[0]);
			int y = getCy(tab[1]);

			if (i != traitPosName) {
				Font taFont = new Font("SansSerif", Font.PLAIN, 10);
				g.setFont(taFont);
			}

			g.setColor(nameColor);
			g.drawString(temp, x, y);

			Font defaut = new Font("SansSerif", Font.PLAIN, 12);
			g.setFont(defaut);

		}

	}

	public void addWaterMediated(Path p) {
		currentWaterMediated.add(p);
		traitPath(p);
		repaint();
	}

	public void deleteWaterMediated(Path p) {
		for (int i = 0; i < currentWaterMediated.size(); i++) {
			Path pp = (Path) currentWaterMediated.get(i);
			if (pp.isEqual(p)) {
				currentWaterMediated.remove(i);
				repaint();
				return;
			}
		}
	}

	public void paintName2(Graphics g) {

		for (int i = 0; i < resExplicit.size(); i++) {

			Explicit e = (Explicit) resExplicit.get(i);
			LinkedList tmp = e.getResidues();
			Residue r = (Residue) tmp.get(0);

			if (r.getType().equals("ion"))
				return;

			Chain c = r.getParent();

			String temp = r.getName() + " " + r.getId() + "(" + c.getName()
					+ ")";
			double[] tab = (double[]) posName2.get(i);

			int x = getCx(tab[0]);
			int y = getCy(tab[1]);

			if (i != traitPosName2) {
				Font taFont = new Font("SansSerif", Font.PLAIN, 10);
				g.setFont(taFont);
			}

			g.setColor(name2Color);
			g.drawString(temp, x, y);

			Font defaut = new Font("SansSerif", Font.PLAIN, 12);
			g.setFont(defaut);

		}

	}

	public void setCutoffC(double cut) {

		if (lockedMode == 0) {
			ligand = new Ligand(ligand.getParent(), cut, cutoffHB);

			atoms = ligand.getAtoms();
			bonds = ligand.getBonds();
			Hbonds = ligand.getHbond();
			resExplicit = ligand.getResExpli();
			explicit = ligand.getExplicit();
			residues = ligand.getResidues();
			contactList = ligand.getContacts();
			rings = ligand.getRing();
			posName = ligand.returnPoseName();
			posName2 = ligand.returnPoseName2();
			init();
			setScale();
		}

		cutoffC = cut;

	}

	public void paintContact(Graphics g, Contact con) {

		double dd = con.getDistance() / 2;
		if (dd <= cutoffC) {

			Atom atom = con.getTrait();
			Residue r = con.getResidue();
			String tmp = "";
			tmp = tmp + r.getId();
			tmp = tmp + "(" + r.getParent().getName() + ")";

			// atom.setXg(getCx(atom.getFx()));
			// atom.setYg(getCy(atom.getFy()));

			int x = getCx(atom.getFx());
			int y = getCy(atom.getFy());

			g.setColor(getResColor(colorMode, r.getName(), r));
			g.fillOval(x - 5, y - 15, 40, 30);
			if (con.getSelected()) {
				g.setColor(Color.RED);
			} else {
				g.setColor(Color.BLACK);
			}
			g.drawOval(x - 5, y - 15, 40, 30);
			g.setColor(Color.BLACK);

			Font taFont = new Font("SansSerif", Font.PLAIN, 10);
			g.setFont(taFont);

			g.setColor(Color.BLACK);
			g.drawString(r.getName(), x, y);
			g.drawString(tmp, x, y + 10);

			Font defaut = new Font("SansSerif", Font.PLAIN, 12);
			g.setFont(defaut);

		}

	}

	public Color getResColor(int mode, String r, Residue res) {

		if (mode == 1) {

			if ((r.equals("ARG")) || (r.equals("ASN")) || (r.equals("ASP"))
					|| (r.equals("GLN")) || (r.equals("GLU"))
					|| (r.equals("HIS")) || (r.equals("LYS"))
					|| (r.equals("SER")) || (r.equals("THR")))
				return Color.PINK;

			if ((r.equals("ALA")) || (r.equals("CYS")) || (r.equals("ILE"))
					|| (r.equals("LEU")) || (r.equals("MET"))
					|| (r.equals("PHE")) || (r.equals("PRO"))
					|| (r.equals("TRP")) || (r.equals("TYR"))
					|| (r.equals("VAL")))
				return Color.GREEN;

			return Color.LIGHT_GRAY;

		}

		if (mode == 2) {

			if ((r.equals("ASP")) || (r.equals("GLU")))
				return (new Color(230, 10, 10));
			if ((r.equals("CYS")) || (r.equals("MET")))
				return (new Color(230, 230, 0));
			if ((r.equals("LYS")) || (r.equals("ARG")))
				return (new Color(20, 90, 255));
			if ((r.equals("SER")) || (r.equals("THR")))
				return (new Color(250, 150, 0));
			if ((r.equals("PHE")) || (r.equals("TYR")))
				return (new Color(50, 50, 170));
			if ((r.equals("ASN")) || (r.equals("GLN")))
				return (new Color(0, 220, 220));
			if (r.equals("GLY"))
				return (new Color(235, 235, 235));
			if ((r.equals("LEU")) || (r.equals("VAL")) || (r.equals("ILE")))
				return (new Color(15, 130, 15));
			if (r.equals("ALA"))
				return (new Color(200, 200, 200));
			if (r.equals("TRP"))
				return (new Color(180, 90, 180));
			if (r.equals("HIS"))
				return (new Color(130, 130, 210));
			if (r.equals("PRO"))
				return (new Color(220, 150, 130));

			if (r.equals("A"))
				return (new Color(160, 160, 255));
			if (r.equals("C"))
				return (new Color(255, 140, 75));
			if (r.equals("G"))
				return (new Color(255, 112, 112));
			if (r.equals("T"))
				return (new Color(160, 255, 160));
			if (r.equals("U"))
				return (new Color(184, 184, 184));

			return Color.WHITE;

		}

		if (mode == 3) {

			/* non-polar */
			if ((r.equals("ALA")) || (r.equals("VAL")) || (r.equals("LEU"))
					|| (r.equals("PRO")) || (r.equals("MET"))
					|| (r.equals("PHE")) || (r.equals("ILE"))
					|| (r.equals("TRP")) || (r.equals("GLY")))
				return (new Color(235, 235, 235));

			/* polar (-) */
			if ((r.equals("ASP")) || (r.equals("GLU")))
				return (new Color(230, 10, 10));

			/* polar (+) */
			if ((r.equals("LYS")) || (r.equals("ARG")) || (r.equals("HIS")))
				return (new Color(50, 50, 170));

			/* polar (uncharged) */
			if ((r.equals("SER")) || (r.equals("THR")) || (r.equals("ASN"))
					|| (r.equals("GLN")) || (r.equals("CYS"))
					|| (r.equals("TYR")))
				return (new Color(180, 90, 180));

			return Color.WHITE;

		}

		if (mode == 4) {

			Chain cc = res.getParent();
			Molecule mol = cc.getParent();
			LinkedList struct = mol.getStruct();

			for (int i = 0; i < struct.size(); i++) {
				Structure s = (Structure) struct.get(i);
				char type = s.getType();
				if (s.containRes(res)) {
					if (type == 'H')
						return (new Color(240, 0, 128));
					if (type == 'E')
						return (new Color(255, 255, 0));

				}

				return Color.WHITE;

			}

		}
		if (mode == 5) {
			return resColor;
		}

		return Color.LIGHT_GRAY;

	}

	public void setColorMode(int mode) {
		colorMode = mode;
	}

	public void paintHbond(Graphics g, Hbond h) {

		double dd = h.getDistance();
		if (dd <= cutoffHB) {

			Atom don = h.getDon();
			Atom acc = h.getAcc();

			int x1 = getCx(don.getFx()) + 5;
			int y1 = getCy(don.getFy()) + 10;

			int x2 = getCx(acc.getFx()) + 5;
			int y2 = getCy(acc.getFy()) + 10;

			LinkedList l1 = getTrans(x1, y1, x2, y2);
			LinkedList l2 = getTrans(x2, y2, x1, y1);

			double xx1 = (Double) l1.get(0);
			double yy1 = (Double) l1.get(1);
			double xx2 = (Double) l2.get(0);
			double yy2 = (Double) l2.get(1);

			g.setColor(hbColor);

			g.drawLine((int) xx1, (int) yy1, (int) xx2, (int) yy2);

			if (hbStyle == 2) {
				drawArrow(g, x2, y2, x1, y1);
			}

			/* Display Hbond distances */
			if (displayDist == 1) {
				DecimalFormat df = new DecimalFormat("########.0");
				String distance = df.format(dd);

				int mx = (x1 - 5 + x2 - 5) / 2;
				int my = (y1 - 10 + y2 - 10) / 2;

				g.setColor(Color.white);
				g.fillRect(mx, my, 20, 10);
				g.setColor(Color.white);
				g.drawRect(mx, my, 20, 10);

				g.setColor(hbColor);

				Font taFont = new Font("SansSerif", Font.PLAIN, 10);
				g.setFont(taFont);

				g.drawString(distance, mx - 5, my + 10);

				Font defaut = new Font("SansSerif", Font.PLAIN, 12);
				g.setFont(defaut);
			}

			Atom aRes = h.getAtomRes();

			int x = getCx(aRes.getFx());
			int y = getCy(aRes.getFy());

			if (h.getSelected()) {
				g.setColor(Color.RED);
			} else {
				g.setColor(hbColor);
			}

			int leng;
			if (aRes.getName().length() == 1) {
				leng = 6;
			} else
				leng = 1;
			g.drawString(aRes.getName(), x + leng, y + 10);

			Font taFont = new Font("SansSerif", Font.PLAIN, 10);
			g.setFont(taFont);

			Residue rr = aRes.getParent();
			Chain c = rr.getParent();

			String tmp = rr.getName();
			String tmp2 = rr.getId() + "(" + c.getName() + ")";

			g.drawString(tmp, x + leng, y + 20);
			g.drawString(tmp2, x + leng, y + 30);

			Font defaut = new Font("SansSerif", Font.PLAIN, 12);
			g.setFont(defaut);
		}

	}

	public void paintPath(Graphics g, Path p) {

		Atom don = p.getDep();

		Atom acc = p.getCoord();

		int x1 = getCx(don.getFx()) + 5;
		int y1 = getCy(don.getFy()) + 10;

		int x2 = getCx(acc.getFx()) + 5;
		int y2 = getCy(acc.getFy()) + 10;

		LinkedList l1 = getTrans(x1, y1, x2, y2);
		LinkedList l2 = getTrans(x2, y2, x1, y1);

		double xx1 = (Double) l1.get(0);
		double yy1 = (Double) l1.get(1);
		double xx2 = (Double) l2.get(0);
		double yy2 = (Double) l2.get(1);

		if (p.getSelected()) {
			g.setColor(Color.RED);
		} else {
			g.setColor(wmColor);
		}

		g.drawLine((int) xx1, (int) yy1, (int) xx2, (int) yy2);

		Atom aRes = p.getCoord();

		int x = getCx(aRes.getFx());
		int y = getCy(aRes.getFy());

		Font taFont = new Font("SansSerif", Font.PLAIN, 10);
		g.setFont(taFont);

		String tmp = p.getString2();
		StringTokenizer tok = new StringTokenizer(tmp, ";");
		int c = tok.countTokens();

		for (int i = 1; i <= c; i++) {

			g.drawString(tok.nextToken(), x + 5, y + (i * 10));

		}

		Font defaut = new Font("SansSerif", Font.PLAIN, 12);
		g.setFont(defaut);

	}

	public static void drawArrow(Graphics g, int x1, int y1, int x2, int y2) {

		int xa = 0;
		int ya = 0;
		int xb = 0;
		int yb = 0;

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;
		double r = 10;

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

		double d1 = Math.sqrt((sol1x - x2) * (sol1x - x2) + (sol1y - y2)
				* (sol1y - y2));
		double d2 = Math.sqrt((sol2x - x2) * (sol2x - x2) + (sol2y - y2)
				* (sol2y - y2));

		if (d1 < d2) {
			solx = sol1x;
			soly = sol1y;
		} else {
			solx = sol2x;
			soly = sol2y;
		}

		double sx = solx + (solx - x1);
		double sy = soly + (soly - y1);

		double OAx = sx - solx;
		double OAy = sy - soly;

		double resx = solx + Math.cos(0.69) * OAx - Math.sin(0.69) * OAy;
		double resy = soly + Math.sin(0.69) * OAx + Math.cos(0.69) * OAy;

		double resx2 = solx + Math.cos(-0.69) * OAx - Math.sin(-0.69) * OAy;
		double resy2 = soly + Math.sin(-0.69) * OAx + Math.cos(-0.69) * OAy;

		xa = (int) resx;
		xb = (int) resx2;
		ya = (int) resy;
		yb = (int) resy2;

		int abcisses[] = new int[] { (int) solx, xa, xb };
		int ordonnes[] = new int[] { (int) soly, ya, yb };
		g.setColor(hbColor);
		g.fillPolygon(abcisses, ordonnes, 3);

	}

	public LinkedList getTrans(int x1, int y1, int x2, int y2) {

		double sol1x = 0;
		double sol1y = 0;
		double sol2x = 0;
		double sol2y = 0;

		double a = y2 - y1;
		double b = x1 - x2;
		double c = y1 * x2 - x1 * y2;

		double xo = x1;
		double yo = y1;
		double r = 10;

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

		double d1 = Math.sqrt((sol1x - x2) * (sol1x - x2) + (sol1y - y2)
				* (sol1y - y2));
		double d2 = Math.sqrt((sol2x - x2) * (sol2x - x2) + (sol2y - y2)
				* (sol2y - y2));

		if (d1 < d2) {
			solx = sol1x;
			soly = sol1y;
		} else {
			solx = sol2x;
			soly = sol2y;
		}

		LinkedList res = new LinkedList();
		res.add(solx);
		res.add(soly);

		return res;

	}

	public void paintEx(Graphics g, Explicit e) {

		Graphics2D g2d = (Graphics2D) g;

		LinkedList bondsE = e.getBonds();
		LinkedList atomsE = e.getAtoms();
		LinkedList ringsE = e.getRing();

		/* bonds */
		for (int i = 0; i < bondsE.size(); i++) {
			Bond b = (Bond) bondsE.get(i);
			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();

			int x1 = getCx(a1.getFx());
			int y1 = getCy(a1.getFy());
			int x2 = getCx(a2.getFx());
			int y2 = getCy(a2.getFy());

			g2d.setStroke(new BasicStroke(2.0f));

			if (e.getSelected()) {
				g.setColor(Color.RED);
			} else {
				g.setColor(resExplicitColor);
			}

			if (b.getOrder() == 1) {
				g.drawLine(x1 + 10, y1 + 5, x2 + 10, y2 + 5);
			}
			if (b.getOrder() == 2) {
				LinkedList l1 = getDouble(a1, a2);
				LinkedList l2 = getDouble(a2, a1);

				int X1 = 0;
				int Y1 = 0;
				int X2 = 0;
				int Y2 = 0;
				int X3 = 0;
				int Y3 = 0;
				int X4 = 0;
				int Y4 = 0;

				double xx1 = (Double) l1.get(0);
				double xx2 = (Double) l1.get(2);
				double yy1 = (Double) l1.get(1);
				double yy2 = (Double) l1.get(3);

				if (xx1 != xx2) {
					if (xx1 < xx2) {
						X1 = getCx(xx1);
						Y1 = getCy((Double) l1.get(1));
						X2 = getCx(xx2);
						Y2 = getCy((Double) l1.get(3));
					} else {
						X1 = getCx(xx2);
						Y1 = getCy((Double) l1.get(3));
						X2 = getCx(xx1);
						Y2 = getCy((Double) l1.get(1));
					}
				} else {
					if (yy1 < yy2) {
						X1 = getCx(xx1);
						Y1 = getCy((Double) l1.get(1));
						X2 = getCx(xx2);
						Y2 = getCy((Double) l1.get(3));
					} else {
						X1 = getCx(xx2);
						Y1 = getCy((Double) l1.get(3));
						X2 = getCx(xx1);
						Y2 = getCy((Double) l1.get(1));
					}
				}

				double xxx1 = (Double) l2.get(0);
				double xxx2 = (Double) l2.get(2);
				double yyy1 = (Double) l2.get(1);
				double yyy2 = (Double) l2.get(3);

				if (xxx1 != xxx2) {
					if (xxx1 < xxx2) {
						X3 = getCx(xxx1);
						Y3 = getCy((Double) l2.get(1));
						X4 = getCx(xxx2);
						Y4 = getCy((Double) l2.get(3));
					} else {
						X3 = getCx(xxx2);
						Y3 = getCy((Double) l2.get(3));
						X4 = getCx(xxx1);
						Y4 = getCy((Double) l2.get(1));
					}
				} else {
					if (yyy1 < yyy2) {
						X3 = getCx(xxx1);
						Y3 = getCy((Double) l2.get(1));
						X4 = getCx(xxx2);
						Y4 = getCy((Double) l2.get(3));
					} else {
						X3 = getCx(xxx2);
						Y3 = getCy((Double) l2.get(3));
						X4 = getCx(xxx1);
						Y4 = getCy((Double) l2.get(1));
					}
				}

				g.drawLine(X1 + 10, Y1 + 5, X3 + 10, Y3 + 5);
				g.drawLine(X2 + 10, Y2 + 5, X4 + 10, Y4 + 5);
			}
			if (b.getOrder() == 3) {

				g.drawLine(x1 + 10, y1 + 5, x2 + 10, y2 + 5);

				LinkedList l1 = getDouble(a1, a2);
				LinkedList l2 = getDouble(a2, a1);

				int X1 = 0;
				int Y1 = 0;
				int X2 = 0;
				int Y2 = 0;
				int X3 = 0;
				int Y3 = 0;
				int X4 = 0;
				int Y4 = 0;

				double xx1 = (Double) l1.get(0);
				double xx2 = (Double) l1.get(2);
				double yy1 = (Double) l1.get(1);
				double yy2 = (Double) l1.get(3);

				if (xx1 != xx2) {
					if (xx1 < xx2) {
						X1 = getCx(xx1);
						Y1 = getCy((Double) l1.get(1));
						X2 = getCx(xx2);
						Y2 = getCy((Double) l1.get(3));
					} else {
						X1 = getCx(xx2);
						Y1 = getCy((Double) l1.get(3));
						X2 = getCx(xx1);
						Y2 = getCy((Double) l1.get(1));
					}
				} else {
					if (yy1 < yy2) {
						X1 = getCx(xx1);
						Y1 = getCy((Double) l1.get(1));
						X2 = getCx(xx2);
						Y2 = getCy((Double) l1.get(3));
					} else {
						X1 = getCx(xx2);
						Y1 = getCy((Double) l1.get(3));
						X2 = getCx(xx1);
						Y2 = getCy((Double) l1.get(1));
					}
				}

				double xxx1 = (Double) l2.get(0);
				double xxx2 = (Double) l2.get(2);
				double yyy1 = (Double) l2.get(1);
				double yyy2 = (Double) l2.get(3);

				if (xxx1 != xxx2) {
					if (xxx1 < xxx2) {
						X3 = getCx(xxx1);
						Y3 = getCy((Double) l2.get(1));
						X4 = getCx(xxx2);
						Y4 = getCy((Double) l2.get(3));
					} else {
						X3 = getCx(xxx2);
						Y3 = getCy((Double) l2.get(3));
						X4 = getCx(xxx1);
						Y4 = getCy((Double) l2.get(1));
					}
				} else {
					if (yyy1 < yyy2) {
						X3 = getCx(xxx1);
						Y3 = getCy((Double) l2.get(1));
						X4 = getCx(xxx2);
						Y4 = getCy((Double) l2.get(3));
					} else {
						X3 = getCx(xxx2);
						Y3 = getCy((Double) l2.get(3));
						X4 = getCx(xxx1);
						Y4 = getCy((Double) l2.get(1));
					}
				}

				g.drawLine(X1 + 10, Y1 + 5, X3 + 10, Y3 + 5);
				g.drawLine(X2 + 10, Y2 + 5, X4 + 10, Y4 + 5);

			}

		}

		g2d.setStroke(new BasicStroke(1.0f));

		/* atoms */
		if (labelMode == 0) {
			for (int i = 0; i < atomsE.size(); i++) {

				Atom atom = (Atom) atomsE.get(i);

				// atom.setXg(getCx(atom.getFx()));
				// atom.setYg(getCy(atom.getFy()));

				int x = getCx(atom.getFx());
				int y = getCy(atom.getFy());

				g.setColor(Color.white);
				g.fillRect(x, y, 20, 10);
				g.setColor(Color.white);
				g.drawRect(x, y, 20, 10);

				if (e.getSelected()) {
					g.setColor(Color.RED);
				} else {
					if (atomColorMode == 1) {
						g.setColor(resExplicitColor);
					} else {
						g.setColor(getAtomColor(atom.getElement()));
					}
				}

				int leng;
				if (atom.getName().length() == 1) {
					leng = 6;
				} else
					leng = 1;
				g.drawString(atom.getName(), x + leng, y + 10);

			}

		} else {

			for (int i = 0; i < atomsE.size(); i++) {

				Atom atom = (Atom) atomsE.get(i);

				if (!atom.getElement().equals("C")) {

					// atom.setXg(getCx(atom.getFx()));
					// atom.setYg(getCy(atom.getFy()));

					int x = getCx(atom.getFx());
					int y = getCy(atom.getFy());

					g.setColor(Color.white);
					g.fillRect(x + 5, y, 10, 10);
					g.setColor(Color.white);
					g.drawRect(x + 5, y, 10, 10);

					if (e.getSelected()) {
						g.setColor(Color.RED);
					} else {
						if (atomColorMode == 1) {
							g.setColor(resExplicitColor);
						} else {
							g.setColor(getAtomColor(atom.getElement()));
						}
					}

					g.drawString(atom.getElement(), x + 5, y + 10);

				}
				/* Carbon */
				else {
					// atom.setXg(getCx(atom.getFx()));
					// atom.setYg(getCy(atom.getFy()));

					int x = getCx(atom.getFx());
					int y = getCy(atom.getFy());

					if (e.getSelected()) {
						g.setColor(Color.RED);
					} else {
						if (atomColorMode == 1) {
							g.setColor(resExplicitColor);
						} else {
							g.setColor(getAtomColor(atom.getElement()));
						}
					}

					g.drawOval(x + 7, y + 2, 5, 5);
					g.fillOval(x + 7, y + 2, 5, 5);

				}

			}

		}

		if (e.getSelected()) {
			g.setColor(Color.RED);
		} else {
			g.setColor(resExplicitColor);
		}

		/* rings */

		g2d.setStroke(new BasicStroke(2.0f));

		if (ringsE != null) {

			for (int i = 0; i < ringsE.size(); i++) {
				Ring r = (Ring) ringsE.get(i);
				if (r.isAromatic()) {
					LinkedList att = r.getAtoms();
					int id = (Integer) att.get(0);

					Atom aaa = getAtom(id, atomsE);
					LinkedList lll = getCenter(r, atomsE);
					double xo = (Double) lll.get(0);
					double yo = (Double) lll.get(1);

					int x1 = getCx(xo);
					int y1 = getCy(yo);

					int x2 = getCx(aaa.getFx());
					int y2 = getCy(aaa.getFy());

					double d = Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2)
							* (y1 - y2));
					int dd = (int) Math.rint(d);

					g.drawOval(x1 + 10 - (dd / 2), y1 + 5 - (dd / 2), dd, dd);

				}
			}
		}

		g2d.setStroke(new BasicStroke(1.0f));

	}

	public LinkedList getCenter(Ring r, LinkedList atoms) {
		LinkedList l = r.getAtoms();
		LinkedList tmp = new LinkedList();
		for (int i = 0; i < l.size(); i++) {
			int id = (Integer) l.get(i);
			Atom aaa = getAtom(id, atoms);
			tmp.add(aaa);
		}
		return barycenter(tmp);
	}

	public Atom getAtom(int id, LinkedList l) {
		for (int i = 0; i < l.size(); i++) {
			Atom at = (Atom) l.get(i);
			if (at.getId() == id)
				return at;
		}
		return null;
	}

	public void setScale() {

		viewSize = new double[4];
		view = new Atom[4];

		viewSize[0] = 9.9e19;
		viewSize[1] = -9.9e19;
		viewSize[2] = 9.9e19;
		viewSize[3] = -9.9e19;

		for (int i = 0; i < atoms.size(); i++) {
			Atom atom = (Atom) atoms.get(i);
			atom.getSize(viewSize, view);
		}

		for (int j = 0; j < resExplicit.size(); j++) {
			Explicit e = (Explicit) resExplicit.get(j);

			LinkedList atomsE = e.getAtoms();
			for (int i = 0; i < atomsE.size(); i++) {
				Atom atom = (Atom) atomsE.get(i);
				atom.getSize(viewSize, view);
			}
		}

		for (int i = 0; i < Hbonds.size(); i++) {
			Hbond h = (Hbond) Hbonds.get(i);
			Atom atom = h.getAtomRes();
			atom.getSize(viewSize, view);
		}

		for (int i = 0; i < contactList.size(); i++) {
			Contact con = (Contact) contactList.get(i);
			Atom atom = con.getTrait();
			atom.getSize(viewSize, view);
		}

		for (int i = 0; i < currentWaterMediated.size(); i++) {
			Path p = (Path) currentWaterMediated.get(i);
			Atom atom = p.getCoord();
			atom.getSize(viewSize, view);
		}

		for (int i = 0; i < posName.size(); i++) {

		}

		Atom att = (Atom) atoms.get(0);

		for (int i = 0; i < posName.size(); i++) {
			double[] tab = (double[]) posName.get(i);
			if (tab != null) {
				Atom at = att.clone();
				at.setFx(tab[0]);
				at.setFy(tab[1]);
				at.getSize(viewSize, view);
			}
		}

		for (int i = 0; i < posName2.size(); i++) {
			double[] tab = (double[]) posName2.get(i);
			if (tab != null) {
				Atom at = att.clone();
				at.setFx(tab[0]);
				at.setFy(tab[1]);
				at.getSize(viewSize, view);
			}
		}

		viewSize[0] = viewSize[0]  -0.2;
		viewSize[1] = viewSize[1]  ;
		viewSize[2] = viewSize[2]  -0.9;
		viewSize[3] = viewSize[3]  ;

	}

	public LinkedList getAll() {

		LinkedList result = new LinkedList();

		for (int i = 0; i < atoms.size(); i++) {
			Atom atom = (Atom) atoms.get(i);
			result.add(atom);
		}

		for (int j = 0; j < resExplicit.size(); j++) {
			Explicit e = (Explicit) resExplicit.get(j);

			LinkedList atomsE = e.getAtoms();
			for (int i = 0; i < atomsE.size(); i++) {
				Atom atom = (Atom) atomsE.get(i);
				result.add(atom);
			}
		}

		for (int i = 0; i < Hbonds.size(); i++) {
			Hbond h = (Hbond) Hbonds.get(i);
			Atom atom = h.getAtomRes();
			if (!isInside(atom, result)) {
				result.add(atom);
			}
		}

		for (int i = 0; i < contactList.size(); i++) {
			Contact con = (Contact) contactList.get(i);
			Atom atom = con.getTrait();
			result.add(atom);
		}

		return result;

	}

	public boolean isInside(Atom at, LinkedList l) {

		for (int i = 0; i < l.size(); i++) {
			Atom a = (Atom) l.get(i);
			if (a.getId() == at.getId())
				return true;
		}
		return false;
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

	public void rotation(LinkedList l, double x0, double y0, double angle) {

		angle = Math.PI * angle / 180;

		for (int i = 0; i < l.size(); i++) {
			Atom at = (Atom) l.get(i);

			double OAx = at.getFx() - x0;
			double OAy = at.getFy() - y0;

			double resx = x0 + Math.cos(angle) * OAx - Math.sin(angle) * OAy;
			double resy = y0 + Math.sin(angle) * OAx + Math.cos(angle) * OAy;

			at.setFx(resx);
			at.setFy(resy);

		}

		for (int i = 0; i < posName.size(); i++) {
			double[] tab = (double[]) posName.get(i);
			if (tab != null) {
				double OAx = tab[0] - x0;
				double OAy = tab[1] - y0;

				double resx = x0 + Math.cos(angle) * OAx - Math.sin(angle)
						* OAy;
				double resy = y0 + Math.sin(angle) * OAx + Math.cos(angle)
						* OAy;

				tab[0] = resx;
				tab[1] = resy;
			}

		}

		for (int i = 0; i < posName2.size(); i++) {
			double[] tab = (double[]) posName2.get(i);
			if (tab != null) {
				double OAx = tab[0] - x0;
				double OAy = tab[1] - y0;

				double resx = x0 + Math.cos(angle) * OAx - Math.sin(angle)
						* OAy;
				double resy = y0 + Math.sin(angle) * OAx + Math.cos(angle)
						* OAy;

				tab[0] = resx;
				tab[1] = resy;
			}

		}

		for (int i = 0; i < rings.size(); i++) {
			Ring r = (Ring) rings.get(i);
			if (r.isAromatic()) {

				double xo = r.getCenterX();
				double yo = r.getCenterY();

				double OAx = xo - x0;
				double OAy = yo - y0;

				double resx = x0 + Math.cos(angle) * OAx - Math.sin(angle)
						* OAy;
				double resy = y0 + Math.sin(angle) * OAx + Math.cos(angle)
						* OAy;

				r.setCenterX(resx);
				r.setCenterY(resy);

			}
		}

		for (int i = 0; i < currentWaterMediated.size(); i++) {
			Path p = (Path) currentWaterMediated.get(i);
			Atom at = p.getCoord();

			double xo = at.getFx();
			double yo = at.getFy();

			double OAx = xo - x0;
			double OAy = yo - y0;

			double resx = x0 + Math.cos(angle) * OAx - Math.sin(angle) * OAy;
			double resy = y0 + Math.sin(angle) * OAx + Math.cos(angle) * OAy;

			at.setFx(resx);
			at.setFy(resy);

		}

	}

	public void rotation() {

		LinkedList all = getAll();

		LinkedList tmp = barycenter(all);
		rotation(all, (Double) tmp.get(0), (Double) tmp.get(1), -45);
		setScale();
		drawScale3();
		repaint();

	}

	private void drawScale() {

		double[] mapping = new double[4];

		double dx = viewSize[1] - viewSize[0];
		double dy = viewSize[3] - viewSize[2];
		double scale;
		if (dx > dy)
			scale = dx;
		else
			scale = dy;
		if (scale < 0.0001)
			scale = 1.0;

		mapping[1] = viewSize[0];
		mapping[2] = viewSize[2];
		mapping[0] = (size().width - 50) / scale;

		for (int i = 0; i < residues.size(); i++) {
			Residue res = (Residue) residues.get(i);

			res.setMapping(mapping);
		}

		return;

	}

	private void drawScale2() {

		double dx = (size().width) / (viewSize[1] - viewSize[0]);
		double dy = (size().height) / (viewSize[3] - viewSize[2]);

		scale = Math.max(dx, dy) * 0.4;

	}

	private void drawScale3() {

		double ddx = Math.sqrt((view[1].getFx() - view[0].getFx())
				* (view[1].getFx() - view[0].getFx())
				+ (view[1].getFy() - view[0].getFy())
				* (view[1].getFy() - view[0].getFy()));

		double ddy = Math.sqrt((view[3].getFx() - view[2].getFx())
				* (view[3].getFx() - view[2].getFx())
				+ (view[3].getFy() - view[2].getFy())
				* (view[3].getFy() - view[2].getFy()));

		double dx = ddx / (size().width - 100);
		double dy = ddy / (size().height - 100);

		scale = 1 / Math.max(dx, dy);

	}

	int getCx(double x) {

		return (int) Math.rint((x - viewSize[0]) * scale) + 10;
	}

	double getFFx(int x) {
		return (((x - 10) / scale) + viewSize[0]);

	}

	int getCy(double y) {

		return (size().height) - ((int) Math.rint((y - viewSize[2]) * scale) + 20);
		//return (size().height) - ((int) Math.rint((y - viewSize[2]) * (scale-20)) );
	}

	double getFFy(int y) {

		return (((size().height - y) - 20) / scale) + viewSize[2];
	}

	public void drawCircle(Graphics cg, int xCenter, int yCenter, int r) {
		cg.drawOval(xCenter - r, yCenter - r, 2 * r, 2 * r);
	}

	public void mouseClicked(MouseEvent e) {
	}

	public void mouseDragged(MouseEvent e) {

		if (lockedMode == 1) {
			if (mode == 1) {

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					Atom at = traitContact.getTrait();

					at.setFx(getFFx(tmpX));
					
					at.setFy(getFFy(tmpY));
					
					repaint();
				}

			}

			if (mode == 2) {

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					Atom at = traitHbond.getAtomRes();

					at.setFx(getFFx(tmpX));

					at.setFy(getFFy(tmpY));
					repaint();
				}

			}

			if (mode == 3) {

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double dx = getFFx(tmpX) - atomExplicit.getFx();
					double dy = getFFy(tmpY) - atomExplicit.getFy();

					LinkedList l = traitExplicit.getAtoms();

					for (int i = 0; i < l.size(); i++) {
						Atom at = (Atom) l.get(i);

						at.setFx(at.getFx() + dx);
						at.setFy(at.getFy() + dy);
					}

				}

				repaint();
			}

			if (mode == 4) {

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double[] tab = (double[]) posName.get(traitPosName);

					tab[0] = (getFFx(tmpX));
					tab[1] = (getFFy(tmpY));
				}

				repaint();
			}

			if (mode == 5) {

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double[] tab = (double[]) posName2.get(traitPosName2);

					tab[0] = (getFFx(tmpX));
					tab[1] = (getFFy(tmpY));
				}

				repaint();
			}

			if (mode == 6) {

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					Atom at = traitPath.getCoord();

					at.setFx(getFFx(tmpX));

					at.setFy(getFFy(tmpY));
				}
				repaint();
			}

			if (mode == 7) {

				int tmpX = e.getX();
				int tmpY = e.getY();
				//System.out.println("X:"+e.getX());
				//System.out.println("Y:"+e.getY());

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double dx = getFFx(tmpX) - traitLigand.getFx();
					double dy = getFFy(tmpY) - traitLigand.getFy();

					for (int i = 0; i < atoms.size(); i++) {
						Atom at = (Atom) atoms.get(i);

						at.setFx(at.getFx() + dx);
						at.setFy(at.getFy() + dy);
					}
					for (int i = 0; i < rings.size(); i++) {
						Ring r = (Ring) rings.get(i);
						if (r.isAromatic()) {

							double xo = r.getCenterX();
							double yo = r.getCenterY();

							r.setCenterX(xo + dx);
							r.setCenterY(yo + dy);

						}
					}

				}

				repaint();
			}


			if (mode == 8) {

				int tmpX = e.getX();
				int tmpY = e.getY();
				//System.out.println("X:"+e.getX());
				//System.out.println("Y:"+e.getY());

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double dx = getFFx(tmpX) - traitLigand.getFx();
					double dy = getFFy(tmpY) - traitLigand.getFy();

					/*for (int i = 0; i < atoms.size(); i++) {
						Atom at = (Atom) atoms.get(i);*/

						selectedAtom.setFx(selectedAtom.getFx() + dx);
						selectedAtom.setFy(selectedAtom.getFy() + dy);
					//}
					/*for (int i = 0; i < rings.size(); i++) {
						Ring r = (Ring) rings.get(i);
						if (r.isAromatic()) {

							double xo = r.getCenterX();
							double yo = r.getCenterY();

							r.setCenterX(xo + dx);
							r.setCenterY(yo + dy);

						}
					}*/

				}

				repaint();
			}

		}
	}

	public void mouseMoved(MouseEvent e) {

	}



 public void mousePressed(MouseEvent e) {
	
	/* Mouse left Button */
	if (SwingUtilities.isLeftMouseButton( e )){
	if(lockedMode==1){
	x=e.getX();
	y=e.getY();
	//System.out.println("move "+x+" "+y);
	Contact c=getContact(x,y);
	if(c!=null){
	c.setSelected(true);
	mode=1;
	traitContact=c;
	repaint();
	return;
	}

	Hbond h=getHbond(x,y);
	if(h!=null){
	h.setSelected(true);
	mode=2;
	traitHbond=h;
	repaint();
	return;

	}


	Explicit ex=getExplicit(x,y);
	if(ex!=null){
	ex.setSelected(true);
	mode=3;
	traitExplicit=ex;
	repaint();
	return;

	}

	int v=getPosName(x,y);
	if(v!=-1){
	mode=4;
	traitPosName=v;
	repaint();
	return;

	}

	
	int u=getPosName2(x,y);
	if(u!=-1){
	mode=5;
	traitPosName2=u;
	repaint();
	return;

	}

	Path p=getPath(x,y);
	if(p!=null){
	mode=6;
	p.setSelected(true);
	traitPath=p;
	repaint();
	return;

	}

	
	Atom aaa=getLigand(x,y);
	
	if(aaa!=null){
	System.out.println("atomselected "+aaa.getName());
	ligand.setSelected(true);
	selectedAtom=aaa;
	mode=7;
	repaint();
	return;

	}

	}
	}

	if (SwingUtilities.isRightMouseButton( e )){

	x=e.getX();
	y=e.getY();

	final Hbond h=getHbond(x,y);
	if(h!=null){
	//h.setSelected(true);
	//mode=2;
	//traitHbond=h;
	
	

		JPopupMenu jpopupmenu = new JPopupMenu();

		JMenuItem jmenuitem1 = new JMenuItem("delete" );
            jmenuitem1.addActionListener(new ActionListener() {

                 public void actionPerformed(ActionEvent actionevent)
                 {
                    deleteHbond(h);
			repaint();
                }

            });
	             jpopupmenu.add(jmenuitem1);
			this.add(jpopupmenu);
			//if( event.isPopupTrigger() )
		//{
			jpopupmenu.show( e.getComponent(),e.getX(), e.getY() );
		//}

	
	return;

	}



	final Contact c=getContact(x,y);
	if(c!=null){
	//h.setSelected(true);
	//mode=2;
	//traitHbond=h;
	
	

		JPopupMenu jpopupmenu = new JPopupMenu();

		JMenuItem jmenuitem1 = new JMenuItem("delete" );
            jmenuitem1.addActionListener(new ActionListener() {

                 public void actionPerformed(ActionEvent actionevent)
                 {
                    deleteContact(c);
			repaint();
                }

            });
	             jpopupmenu.add(jmenuitem1);
			this.add(jpopupmenu);
			//if( event.isPopupTrigger() )
		//{
			jpopupmenu.show( e.getComponent(),e.getX(), e.getY() );
		//}

	
	return;

	}
	



	}

	}

	public void mouseReleased(MouseEvent e) {

		if (lockedMode == 1) {
			if (mode == 1) {

				mode = 0;
				traitContact.setSelected(false);

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					Atom at = traitContact.getTrait();

					at.setFx(getFFx(tmpX));
					at.setFy(getFFy(tmpY));

					traitContact = null;

				}

				repaint();
			}

			if (mode == 2) {

				mode = 0;
				traitHbond.setSelected(false);

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					Atom at = traitHbond.getAtomRes();

					at.setFx(getFFx(tmpX));
					at.setFy(getFFy(tmpY));

					traitHbond = null;

				}

				repaint();
			}

			if (mode == 3) {

				mode = 0;
				traitExplicit.setSelected(false);

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double dx = getFFx(tmpX) - atomExplicit.getFx();
					double dy = getFFy(tmpY) - atomExplicit.getFy();

					LinkedList l = traitExplicit.getAtoms();

					for (int i = 0; i < l.size(); i++) {
						Atom at = (Atom) l.get(i);

						at.setFx(at.getFx() + dx);
						at.setFy(at.getFy() + dy);
					}

					traitExplicit = null;
					atomExplicit = null;

				}

				repaint();
			}

			if (mode == 4) {

				mode = 0;

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double[] tab = (double[]) posName.get(traitPosName);

					tab[0] = (getFFx(tmpX));
					tab[1] = (getFFy(tmpY));

					traitPosName = -1;

				}

				repaint();
			}

			if (mode == 5) {

				mode = 0;

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double[] tab = (double[]) posName2.get(traitPosName2);

					tab[0] = (getFFx(tmpX));
					tab[1] = (getFFy(tmpY));

					traitPosName2 = -1;

				}

				repaint();
			}

			if (mode == 6) {

				mode = 0;
				traitPath.setSelected(false);

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					Atom at = traitPath.getCoord();

					at.setFx(getFFx(tmpX));
					at.setFy(getFFy(tmpY));

					traitPath = null;

				}

				repaint();
			}

			if (mode == 7) {
				mode = 0;

				ligand.setSelected(false);

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double dx = getFFx(tmpX) - traitLigand.getFx();
					double dy = getFFy(tmpY) - traitLigand.getFy();

					for (int i = 0; i < atoms.size(); i++) {
						Atom at = (Atom) atoms.get(i);

						at.setFx(at.getFx() + dx);
						at.setFy(at.getFy() + dy);
					}

					for (int i = 0; i < rings.size(); i++) {
						Ring r = (Ring) rings.get(i);
						if (r.isAromatic()) {

							double xo = r.getCenterX();
							double yo = r.getCenterY();

							r.setCenterX(xo + dx);
							r.setCenterY(yo + dy);

						}
					}

					traitLigand = null;

				}

				repaint();
			}

		}

	}

	public void mouseEntered(MouseEvent e) {

	}

	public void mouseExited(MouseEvent e) {
	}

	Contact getContact(int x0, int y0) {

		for (int i = 0; i < contactList.size(); i++) {
			Contact con = (Contact) contactList.get(i);
			Atom at = con.getTrait();

			int x1 = getCx(at.getFx()) - 5;
			int x2 = x1 + 40;
			int y1 = getCy(at.getFy()) - 15;
			int y2 = y1 + 30;

			if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2))
				return con;

		}

		return null;
	}

	Path getPath(int x0, int y0) {

		for (int i = 0; i < currentWaterMediated.size(); i++) {
			Path p = (Path) currentWaterMediated.get(i);
			Atom at = p.getCoord();

			int x1 = getCx(at.getFx()) - 5;
			int x2 = x1 + 40;
			int y1 = getCy(at.getFy()) - 10;
			int y2 = y1 + 30;

			if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2))
				return p;

		}

		return null;
	}

	int getPosName(int x0, int y0) {

		for (int i = 0; i < posName.size(); i++) {

			double[] tab = (double[]) posName.get(i);

			int x1 = getCx(tab[0]);
			int x2 = x1 + 50;
			int y1 = getCy(tab[1]) - 10;
			int y2 = y1 + 15;

			if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2))
				return i;

		}

		return -1;
	}

	int getPosName2(int x0, int y0) {

		for (int i = 0; i < posName2.size(); i++) {

			double[] tab = (double[]) posName2.get(i);
			if (tab == null)
				return -1;

			int x1 = getCx(tab[0]);
			int x2 = x1 + 50;
			int y1 = getCy(tab[1]) - 10;
			int y2 = y1 + 15;

			if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2))
				return i;

		}

		return -1;
	}

	Hbond getHbond(int x0, int y0) {

		Hbond hr = null;
		for (int i = 0; i < Hbonds.size(); i++) {
			Hbond h = (Hbond) Hbonds.get(i);
			Atom at = h.getAtomRes();

			int x1 = getCx(at.getFx());
			int x2 = x1 + 40;
			int y1 = getCy(at.getFy());
			int y2 = y1 + 40;

			if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2))
				hr = h;

		}

		return hr;
	}

	Explicit getExplicit(int x0, int y0) {

		for (int j = 0; j < resExplicit.size(); j++) {
			Explicit ee = (Explicit) resExplicit.get(j);

			LinkedList atomsE = ee.getAtoms();
			for (int i = 0; i < atomsE.size(); i++) {
				Atom atom = (Atom) atomsE.get(i);

				int x1 = getCx(atom.getFx());
				int x2 = x1 + 50;
				int y1 = getCy(atom.getFy());
				int y2 = y1 + 50;

				if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2)) {
					atomExplicit = atom;
					return ee;
				}

			}
		}

		return null;
	}

	Atom getLigand(int x0, int y0) {

		for (int i = 0; i < atoms.size(); i++) {
			Atom atom = (Atom) atoms.get(i);

			int x1 = getCx(atom.getFx());
			int x2 = x1 + 50;
			int y1 = getCy(atom.getFy());
			int y2 = y1 + 50;

			if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2)) {
				traitLigand = atom;
				return atom;
			}

		}

		return null;
	}

	public void traitPath(Path p) {

		Atom ref = p.getDep();

		Atom res = p.getCoord();

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

		double r = p.getDist() * 2;

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

		res.setFx(solx);
		res.setFy(soly);

		if (conflit(p)) {

			double[] rot = new double[73];

			rot[0] = 5;
			rot[1] = 5;
			rot[2] = 5;
			rot[3] = 5;
			rot[4] = 5;
			rot[5] = 5;
			rot[6] = 5;
			rot[7] = 5;
			rot[8] = 5;
			rot[9] = 5;
			rot[10] = 5;
			rot[11] = 5;
			rot[12] = 5;
			rot[13] = 5;
			rot[14] = 5;
			rot[15] = 5;
			rot[16] = 5;
			rot[17] = 5;
			rot[18] = 5;
			rot[19] = 5;
			rot[20] = 5;
			rot[21] = 5;
			rot[22] = 5;
			rot[23] = 5;
			rot[24] = 5;
			rot[25] = 5;
			rot[26] = 5;
			rot[27] = 5;
			rot[28] = 5;
			rot[29] = 5;
			rot[30] = 5;
			rot[31] = 5;
			rot[32] = 5;
			rot[33] = 5;
			rot[34] = 5;
			rot[35] = 5;
			rot[36] = 5;
			rot[37] = 5;
			rot[38] = 5;
			rot[39] = 5;
			rot[40] = 5;
			rot[41] = 5;
			rot[42] = 5;
			rot[43] = 5;
			rot[44] = 5;
			rot[45] = 5;
			rot[46] = 5;
			rot[47] = 5;
			rot[48] = 5;
			rot[49] = 5;
			rot[50] = 5;
			rot[51] = 5;
			rot[52] = 5;
			rot[53] = 5;
			rot[54] = 5;
			rot[55] = 5;
			rot[56] = 5;
			rot[57] = 5;
			rot[58] = 5;
			rot[59] = 5;
			rot[60] = 5;
			rot[61] = 5;
			rot[62] = 5;
			rot[63] = 5;
			rot[64] = 5;
			rot[65] = 5;
			rot[66] = 5;
			rot[67] = 5;
			rot[68] = 5;
			rot[69] = 5;
			rot[70] = 5;
			rot[71] = 5;
			rot[72] = -360;

			for (int v = 0; v < rot.length; v++) {

				if (conflit(p)) {

					rotation(p, rot[v]);
				} else {

					break;
				}

			}

		}

	}

	public void rotation(Path p, double angle) {

		angle = angle * Math.PI / 180;

		Atom at = p.getCoord();
		Atom pivot = p.getDep();

		double OAx = at.getFx() - pivot.getFx();
		double OAy = at.getFy() - pivot.getFy();

		double resx = pivot.getFx() + Math.cos(angle) * OAx - Math.sin(angle)
				* OAy;
		double resy = pivot.getFy() + Math.sin(angle) * OAx + Math.cos(angle)
				* OAy;

		at.setFx(resx);
		at.setFy(resy);

	}

	public boolean conflit(Path p) {

		Atom dep = p.getDep();
		Atom end = p.getCoord();

		for (int i = 0; i < atoms.size(); i++) {
			Atom at = (Atom) atoms.get(i);
			if (end.distanceF(at) < 2)
				return true;
		}

		for (int i = 0; i < bonds.size(); i++) {
			Bond b = (Bond) bonds.get(i);
			if (isIntersect(dep, end, b))
				return true;
		}

		for (int s = 0; s < resExplicit.size(); s++) {
			Explicit e = (Explicit) resExplicit.get(s);
			LinkedList atomsE = e.getAtoms();
			LinkedList bondsE = e.getBonds();

			for (int i = 0; i < atomsE.size(); i++) {
				Atom at = (Atom) atomsE.get(i);
				if (end.distanceF(at) < 2)
					return true;
			}

			for (int i = 0; i < bondsE.size(); i++) {
				Bond b = (Bond) bondsE.get(i);
				if (isIntersect(dep, end, b))
					return true;
			}

		}

		for (int s = 0; s < Hbonds.size(); s++) {
			Hbond hb = (Hbond) Hbonds.get(s);
			Atom at = hb.getAtomRes();
			if (end.distanceF(at) < 2)
				return true;
			if (isIntersect(dep, end, hb))
				return true;
		}

		for (int s = 0; s < contactList.size(); s++) {
			Contact con = (Contact) contactList.get(s);
			Atom at = con.getTrait();
			if (end.distanceF(at) < 2)
				return true;
		}

		for (int s = 0; s < currentWaterMediated.size(); s++) {
			Path pp = (Path) currentWaterMediated.get(s);
			if (p != pp) {
				Atom at = pp.getCoord();
				if (end.distanceF(at) < 2)
					return true;
				Bond b = new Bond(pp.getDep(), pp.getCoord(), 1);
				if (isIntersect(dep, end, b))
					return true;
			}
		}

		return false;
	}

	public boolean isIntersect(Atom aa1, Atom aa2, Bond bond2) {

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

	public boolean isIntersect(Atom bb1, Atom bb2, Hbond bond1) {

		Atom aa1 = bond1.getDon();
		Atom aa2 = bond1.getAcc();

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
