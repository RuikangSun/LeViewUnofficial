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
 This class defines the panel for the metal frame.
 */

import java.awt.Color;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.*;

import java.awt.BasicStroke;

import javax.swing.JPanel;

import java.io.*;
import java.util.*;

public class MyPanelMetal extends JPanel implements MouseListener,
		MouseMotionListener {

	Metal metal;
	double[] viewSize;
	LinkedList residues;
	LinkedList explicit;
	LinkedList explicit2;
	Atom[] view;
	LinkedList ligands;

	LinkedList resExplicit;
	LinkedList water;
	LinkedList contactList;
	Atom mainAtom;

	LinkedList posName;
	LinkedList posName2;

	MyFrameMetal frame;
	TopMenuMetal panel;
	static int x, y;

	Color brown;
	double cutoffC;

	Color metalColor;
	Color resExplicitColor;
	Color ligandsColor;
	Color watersColor;
	Color resColor;
	Color nameColor;
	Color name2Color;

	double scale;

	double valueC;

	int colorMode;
	int mode;
	Contact traitContact;
	Explicit traitExplicit;
	LigandLight traitLigand;
	Atom atomExplicit;
	Atom atomLigand;
	int traitPosName;
	int traitPosName2;
	Atom traitAtom;
	Atom traitWater;

	int lockedMode;

	int labelMode;
	int atomColorMode;

	public MyPanelMetal() {
	}

	public MyPanelMetal(MyFrameMetal f, TopMenuMetal t) {
		frame = f;
		panel = t;
		addMouseListener(this);
		addMouseMotionListener(this);
		mode = 0;
		lockedMode = 0;
	}

	public JPanel affiche(Metal metal) {
		this.metal = metal;

		mainAtom = metal.getAtom();
		resExplicit = metal.getResExpli();
		explicit = metal.getExplicit();
		explicit2 = metal.getExplicit2();
		ligands = metal.getLightLigands();
		contactList = metal.getContacts();
		water = metal.getWaterList();
		brown = new Color(139, 71, 38);
		posName2 = metal.returnPoseName2();
		posName = new LinkedList();

		for (int j = 0; j < ligands.size(); j++) {

			LigandLight l = (LigandLight) ligands.get(j);

			LinkedList rr = l.getResidues();

			for (int i = 0; i < rr.size(); i++) {

				Residue r = (Residue) rr.get(i);
				Chain c = r.getParent();

				String temp = r.getName() + " " + r.getId() + "(" + c.getName()
						+ ")";
				LinkedList tmp = l.returnPoseName();
				for (int s = 0; s < tmp.size(); s++) {
					posName.add((double[]) tmp.get(s));
				}

			}

		}
		init();

		cutoffC = 4.0;
		traitPosName = -1;
		traitPosName2 = -1;
		colorMode = 1;
		valueC = 4.0;

		labelMode = 1;
		atomColorMode = 0;

		metalColor = Color.GREEN;
		resExplicitColor = Color.BLACK;
		ligandsColor = Color.GRAY;
		watersColor = brown;
		resColor = Color.WHITE;

		nameColor = Color.GRAY;
		name2Color = Color.BLACK;

		setScale();
		this.setBounds(10, 10, 300, 300);
		this.setBackground(Color.WHITE);
		return this;
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		setScale();
		drawScale3();
		paintAtoms(g);
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
		mainAtom = mainAtom.getCopy();

		/* copying resExplicit */

		temp = new LinkedList();
		for (int i = 0; i < resExplicit.size(); i++) {
			Explicit e = (Explicit) resExplicit.get(i);
			Explicit e2 = e.getCopy();
			e2.setAtLigand(mainAtom);
			if (!isInsideEx(e2, temp)) {
				temp.add(e2);
			}
		}
		resExplicit = temp;

		temp = new LinkedList();

		/* copying liagnds */

		for (int i = 0; i < ligands.size(); i++) {
			LigandLight e = (LigandLight) ligands.get(i);
			LigandLight e2 = e.getCopy();
			temp.add(e2);
		}
		ligands = temp;
		temp = new LinkedList();

		for (int i = 0; i < explicit.size(); i++) {
			Bond b = (Bond) explicit.get(i);

			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();
			if (a1.getId() == mainAtom.getId()) {
				Bond b1 = new Bond(mainAtom, getCorres2(a2), b.getOrder());
				temp.add(b1);
			} else {
				Bond b1 = new Bond(getCorres2(a1), mainAtom, b.getOrder());
				temp.add(b1);
			}
		}
		explicit = temp;

		temp = new LinkedList();
		for (int i = 0; i < water.size(); i++) {
			Atom a = (Atom) water.get(i);
			temp.add(a.getCopy());
		}
		water = temp;

		temp = new LinkedList();
		for (int i = 0; i < explicit2.size(); i++) {
			Bond b = (Bond) explicit2.get(i);

			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();
			if (a1.getId() == mainAtom.getId()) {
				Bond b1 = new Bond(mainAtom, getCorres3(a2), b.getOrder());
				temp.add(b1);
			} else {
				Bond b1 = new Bond(getCorres3(a1), mainAtom, b.getOrder());
				temp.add(b1);
			}
		}
		explicit2 = temp;

		temp = new LinkedList();
		for (int i = 0; i < contactList.size(); i++) {
			Contact c = (Contact) contactList.get(i);
			Contact c2 = new Contact(c.getResidue(), c.getTrait().getCopy(), c
					.getRef().getCopy(), c.getDistance() / 2);
			temp.add(c2);
		}
		contactList = temp;

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
		for (int i = 0; i < ligands.size(); i++) {
			LigandLight e = (LigandLight) ligands.get(i);
			LinkedList l = e.getAtoms();
			for (int j = 0; j < l.size(); j++) {
				Atom a1 = (Atom) l.get(j);
				if (at.getId() == a1.getId())
					return a1;
			}
		}
		return null;
	}

	public Atom getCorres3(Atom at) {

		for (int j = 0; j < water.size(); j++) {
			Atom a1 = (Atom) water.get(j);
			if (at.getId() == a1.getId())
				return a1;
		}

		return null;
	}

	public void setCutoffC(double cut) {

		if (lockedMode == 0) {
			metal = new Metal(metal.getParent(), cut);

			resExplicit = metal.getResExpli();
			explicit = metal.getExplicit();
			explicit2 = metal.getExplicit2();
			ligands = metal.getLightLigands();
			contactList = metal.getContacts();
			water = metal.getWaterList();
			mainAtom=metal.getAtom();

			posName2 = metal.returnPoseName2();
			posName = new LinkedList();

			for (int j = 0; j < ligands.size(); j++) {

				LigandLight l = (LigandLight) ligands.get(j);

				LinkedList rr = l.getResidues();

				for (int i = 0; i < rr.size(); i++) {

					Residue r = (Residue) rr.get(i);
					Chain c = r.getParent();

					String temp = r.getName() + " " + r.getId() + "("
							+ c.getName() + ")";
					LinkedList tmp = l.returnPoseName();
					for (int s = 0; s < tmp.size(); s++) {
						posName.add((double[]) tmp.get(s));
					}

				}

			}

			init();
			setScale();
		}

		cutoffC = cut;

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

	public void paintAtoms(Graphics g) {

		/* bonds */
		for (int i = 0; i < explicit.size(); i++) {
			Bond b = (Bond) explicit.get(i);
			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();

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

		/* H2O */
		for (int i = 0; i < explicit2.size(); i++) {
			Bond b = (Bond) explicit2.get(i);
			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();

			int x1 = getCx(a1.getFx());
			int y1 = getCy(a1.getFy());
			int x2 = getCx(a2.getFx());
			int y2 = getCy(a2.getFy());

			g.setColor(watersColor);

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
		}

		int x = getCx(mainAtom.getFx());
		int y = getCy(mainAtom.getFy());

		g.setColor(Color.white);
		g.fillRect(x, y, 20, 10);
		g.setColor(Color.white);
		g.drawRect(x, y, 20, 10);

		if (mainAtom.getSelected()) {
			g.setColor(Color.RED);
		} else {
			g.setColor(metalColor);
		}
		int leng;
		if (mainAtom.getName().length() == 1) {
			leng = 6;
		} else
			leng = 1;
		g.drawString(mainAtom.getName(), x + leng, y + 10);

		for (int i = 0; i < contactList.size(); i++) {

			Contact con = (Contact) contactList.get(i);

			paintContact(g, con);
		}

		for (int i = 0; i < ligands.size(); i++) {
			LigandLight lig = (LigandLight) ligands.get(i);
			paintLig(g, lig);
		}

		for (int i = 0; i < resExplicit.size(); i++) {
			Explicit e = (Explicit) resExplicit.get(i);
			paintEx(g, e);

		}

		for (int i = 0; i < water.size(); i++) {
			Atom w = (Atom) water.get(i);

			int xx = getCx(w.getFx());
			int yy = getCy(w.getFy());

			if (w.getSelected()) {
				g.setColor(Color.RED);
			} else {
				g.setColor(watersColor);
			}

			Font taFont = new Font("SansSerif", Font.PLAIN, 10);
			g.setFont(taFont);

			String tmp = "H2O " + w.getId();
			g.drawString(tmp, xx + 20, yy + 10);

			Font defaut = new Font("SansSerif", Font.PLAIN, 12);
			g.setFont(defaut);

		}

		paintName(g);

		paintName2(g);

	}

	public void paintContact(Graphics g, Contact con) {

		double dd = con.getDistance() / 2;
		if (dd <= cutoffC) {

			Atom atom = con.getTrait();
			Residue r = con.getResidue();
			String tmp = "";
			tmp = tmp + r.getId();
			tmp = tmp + "(" + r.getParent().getName() + ")";

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

			Font taFont = new Font("SansSerif", Font.PLAIN, 10);
			g.setFont(taFont);

			g.setColor(Color.BLACK);
			g.drawString(r.getName(), x, y);
			g.drawString(tmp, x, y + 10);

			Font defaut = new Font("SansSerif", Font.PLAIN, 12);
			g.setFont(defaut);
		}

	}

	public void paintLig(Graphics g, LigandLight lig) {

		Graphics2D g2d = (Graphics2D) g;
		g2d.setStroke(new BasicStroke(2.0f));

		LinkedList bonds = lig.getBonds();
		LinkedList atoms = lig.getAtoms();

		/* bonds */
		for (int i = 0; i < bonds.size(); i++) {
			Bond b = (Bond) bonds.get(i);
			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();

			int x1 = getCx(a1.getFx());
			int y1 = getCy(a1.getFy());
			int x2 = getCx(a2.getFx());
			int y2 = getCy(a2.getFy());

			if (lig.getSelected()) {
				g.setColor(Color.RED);
			} else {
				g.setColor(ligandsColor);
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

				if (lig.getSelected()) {
					g.setColor(Color.RED);
				} else {
					if (atomColorMode == 1) {
						g.setColor(ligandsColor);
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

					int x = getCx(atom.getFx());
					int y = getCy(atom.getFy());

					g.setColor(Color.white);
					g.fillRect(x + 5, y, 10, 10);
					g.setColor(Color.white);
					g.drawRect(x + 5, y, 10, 10);

					if (lig.getSelected()) {
						g.setColor(Color.RED);
					} else {
						if (atomColorMode == 1) {
							g.setColor(ligandsColor);
						} else {
							g.setColor(getAtomColor(atom.getElement()));
						}
					}

					g.drawString(atom.getElement(), x + 5, y + 10);

				}
				/* Carbon */
				else {

					int x = getCx(atom.getFx());
					int y = getCy(atom.getFy());

					if (lig.getSelected()) {
						g.setColor(Color.RED);
					} else {
						if (atomColorMode == 1) {
							g.setColor(ligandsColor);
						} else {
							g.setColor(getAtomColor(atom.getElement()));
						}
					}

					g.drawOval(x + 7, y + 2, 5, 5);
					g.fillOval(x + 7, y + 2, 5, 5);

				}

			}

			if (lig.getSelected()) {
				g.setColor(Color.RED);
			} else {
				g.setColor(ligandsColor);
			}

			g2d.setStroke(new BasicStroke(2.0f));

			/* rings */
			LinkedList rings = lig.getRing();
			if (rings != null) {

				for (int i = 0; i < rings.size(); i++) {
					Ring r = (Ring) rings.get(i);
					if (r.isAromatic()) {
						LinkedList att = r.getAtoms();
						int id = (Integer) att.get(0);

						Atom aaa = getAtom(id, atoms);
						LinkedList lll = getCenter(r, atoms);
						double xo = (Double) lll.get(0);
						double yo = (Double) lll.get(1);

						int x1 = getCx(xo);
						int y1 = getCy(yo);

						int x2 = getCx(aaa.getFx());
						int y2 = getCy(aaa.getFy());

						double d = Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2)
								* (y1 - y2));
						int dd = (int) Math.rint(d);

						g.drawOval(x1 + 10 - (dd / 2), y1 + 5 - (dd / 2), dd,
								dd);

					}
				}
			}

		}

		g2d.setStroke(new BasicStroke(1.0f));

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

	public void setMetalColor(Color c) {
		metalColor = c;
	}

	public void setResExplicitColor(Color c) {
		resExplicitColor = c;
	}

	public void setLigandsColor(Color c) {
		ligandsColor = c;
	}

	public void setWatersColor(Color c) {
		watersColor = c;
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

	public int getLockedMode() {
		return lockedMode;
	}

	public double getValueC() {
		return valueC;
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

	public Color getResColor(String r) {

		if ((r.equals("ARG")) || (r.equals("ASN")) || (r.equals("ASP"))
				|| (r.equals("GLN")) || (r.equals("GLU")) || (r.equals("HIS"))
				|| (r.equals("LYS")) || (r.equals("SER")) || (r.equals("THR")))
			return Color.PINK;

		if ((r.equals("ALA")) || (r.equals("CYS")) || (r.equals("ILE"))
				|| (r.equals("LEU")) || (r.equals("MET")) || (r.equals("PHE"))
				|| (r.equals("PRO")) || (r.equals("TRP")) || (r.equals("TYR"))
				|| (r.equals("VAL")))
			return Color.GREEN;

		return Color.LIGHT_GRAY;
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

	public void paintName(Graphics g) {

		for (int j = 0; j < ligands.size(); j++) {

			LigandLight l = (LigandLight) ligands.get(j);

			LinkedList rr = l.getResidues();

			for (int i = 0; i < rr.size(); i++) {

				Residue r = (Residue) rr.get(i);
				Chain c = r.getParent();

				String temp = r.getName() + " " + r.getId() + "(" + c.getName()
						+ ")";
				LinkedList posName = l.returnPoseName();
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

	}

	public void paintEx(Graphics g, Explicit e) {

		LinkedList bondsE = e.getBonds();
		LinkedList atomsE = e.getAtoms();
		LinkedList ringsE = e.getRing();

		Graphics2D g2d = (Graphics2D) g;
		g2d.setStroke(new BasicStroke(2.0f));

		/* bonds */
		for (int i = 0; i < bondsE.size(); i++) {
			Bond b = (Bond) bondsE.get(i);
			Atom a1 = b.getFirst();
			Atom a2 = b.getSecond();

			int x1 = getCx(a1.getFx());
			int y1 = getCy(a1.getFy());
			int x2 = getCx(a2.getFx());
			int y2 = getCy(a2.getFy());

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

	public void setScale() {

		viewSize = new double[4];
		view = new Atom[4];

		viewSize[0] = 9.9e19;
		viewSize[1] = -9.9e19;
		viewSize[2] = 9.9e19;
		viewSize[3] = -9.9e19;

		mainAtom.getSize(viewSize, view);

		for (int j = 0; j < ligands.size(); j++) {
			LigandLight lig = (LigandLight) ligands.get(j);

			LinkedList atomsE = lig.getAtoms();
			for (int i = 0; i < atomsE.size(); i++) {
				Atom ato = (Atom) atomsE.get(i);
				ato.getSize(viewSize, view);
			}
		}

		for (int j = 0; j < resExplicit.size(); j++) {
			Explicit e = (Explicit) resExplicit.get(j);

			LinkedList atomsE = e.getAtoms();
			for (int i = 0; i < atomsE.size(); i++) {
				Atom ato = (Atom) atomsE.get(i);
				ato.getSize(viewSize, view);
			}
		}

		for (int j = 0; j < water.size(); j++) {
			Atom w = (Atom) water.get(j);

			w.getSize(viewSize, view);
		}

		for (int i = 0; i < contactList.size(); i++) {
			Contact con = (Contact) contactList.get(i);
			Atom aaa = con.getTrait();
			aaa.getSize(viewSize, view);
		}

		for (int i = 0; i < posName.size(); i++) {
			double[] tab = (double[]) posName.get(i);
			if (tab != null) {
				Atom at = mainAtom.clone();
				at.setFx(tab[0]);
				at.setFy(tab[1]);
				at.getSize(viewSize, view);
			}
		}

		for (int i = 0; i < posName2.size(); i++) {
			double[] tab = (double[]) posName2.get(i);
			if (tab != null) {
				Atom at = mainAtom.clone();
				at.setFx(tab[0]);
				at.setFy(tab[1]);
				at.getSize(viewSize, view);
			}
		}

		viewSize[0] = viewSize[0] - 0.2 ;
		viewSize[1] = viewSize[1] ;
		viewSize[2] = viewSize[2] - 0.9 ;
		viewSize[3] = viewSize[3]  ;

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

		double dx = ddx / (size().width - 50);
		double dy = ddy / (size().height - 50);

		scale = 1 / Math.max(dx, dy);

	}

	int getCx(double x) {

		return (int) Math.rint((x - viewSize[0]) * scale) + 10;
	}

	int getCy(double y) {

		return (size().height)
				- ((int) Math.rint((y - viewSize[2]) * scale) + 20);
	}

	public void drawCircle(Graphics cg, int xCenter, int yCenter, int r) {
		cg.drawOval(xCenter - r, yCenter - r, 2 * r, 2 * r);
	}

	public void setLockedMode(int i) {
		lockedMode = i;
		valueC = cutoffC;
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

	double getFFx(int x) {
		return (((x - 10) / scale) + viewSize[0]);

	}

	double getFFy(int y) {

		return (((size().height - y) - 20) / scale) + viewSize[2];
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

					mainAtom.setFx(getFFx(tmpX));

					mainAtom.setFy(getFFy(tmpY));
					repaint();
				}

			}

			if (mode == 2) {

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

					traitWater.setFx(getFFx(tmpX - 20));

					traitWater.setFy(getFFy(tmpY - 10));
					repaint();
				}

			}

			if (mode == 5) {

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double dx = getFFx(tmpX) - atomLigand.getFx();
					double dy = getFFy(tmpY) - atomLigand.getFy();

					LinkedList l = traitLigand.getAtoms();

					for (int i = 0; i < l.size(); i++) {
						Atom at = (Atom) l.get(i);

						at.setFx(at.getFx() + dx);
						at.setFy(at.getFy() + dy);
					}

				}

				repaint();
			}

			if (mode == 6) {

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

			if (mode == 7) {

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

		}

	}

	public void mouseMoved(MouseEvent e) {

	}

	public void mousePressed(MouseEvent e) {

		if (lockedMode == 1) {
			x = e.getX();
			y = e.getY();

			Atom at = getAtom(x, y);
			if (at != null) {
				at.setSelected(true);
				mode = 1;
				traitAtom = at;
				repaint();
				return;
			}

			Contact c = getContact(x, y);
			if (c != null) {
				c.setSelected(true);
				mode = 2;
				traitContact = c;
				repaint();
				return;
			}

			Explicit ex = getExplicit(x, y);
			if (ex != null) {
				ex.setSelected(true);
				mode = 3;
				traitExplicit = ex;
				repaint();
				return;

			}

			Atom w = getWater(x, y);
			if (w != null) {
				w.setSelected(true);
				mode = 4;
				traitWater = w;
				repaint();
				return;

			}

			LigandLight ll = getLigand(x, y);
			if (ll != null) {
				ll.setSelected(true);
				mode = 5;
				traitLigand = ll;
				repaint();
				return;

			}

			int v = getPosName(x, y);
			if (v != -1) {
				mode = 6;
				traitPosName = v;
				repaint();
				return;

			}

			int u = getPosName2(x, y);
			if (u != -1) {
				mode = 7;
				traitPosName2 = u;
				repaint();
				return;

			}

		}

	}

	public void mouseReleased(MouseEvent e) {

		if (lockedMode == 1) {

			if (mode == 1) {

				mode = 0;
				traitAtom.setSelected(false);

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					mainAtom.setFx(getFFx(tmpX));
					mainAtom.setFy(getFFy(tmpY));

					traitAtom = null;

				}

				repaint();
			}

			if (mode == 2) {

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
				traitWater.setSelected(false);

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					traitWater.setFx(getFFx(tmpX - 20));
					traitWater.setFy(getFFy(tmpY - 10));

					traitWater = null;

				}

				repaint();
			}

			if (mode == 5) {

				mode = 0;
				traitLigand.setSelected(false);

				int tmpX = e.getX();
				int tmpY = e.getY();

				if ((tmpX < 10) || (tmpX > size().width - 10)
						|| (tmpY > size().height - 10) || (tmpY < 10)) {

				} else {

					double dx = getFFx(tmpX) - atomLigand.getFx();
					double dy = getFFy(tmpY) - atomLigand.getFy();

					LinkedList l = traitLigand.getAtoms();

					for (int i = 0; i < l.size(); i++) {
						Atom at = (Atom) l.get(i);

						at.setFx(at.getFx() + dx);
						at.setFy(at.getFy() + dy);
					}

					traitLigand = null;
					atomLigand = null;

				}

				repaint();
			}

			if (mode == 6) {

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

			if (mode == 7) {

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

	Atom getWater(int x0, int y0) {

		for (int i = 0; i < water.size(); i++) {
			Atom w = (Atom) water.get(i);

			int x1 = getCx(w.getFx()) + 20;
			int x2 = x1 + 50;
			int y1 = getCy(w.getFy());
			int y2 = y1 + 10;

			if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2))
				return w;

		}
		return null;

	}

	int getPosName(int x0, int y0) {

		for (int i = 0; i < posName.size(); i++) {

			double[] tab = (double[]) posName.get(i);
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

	LigandLight getLigand(int x0, int y0) {

		for (int u = 0; u < ligands.size(); u++) {
			LigandLight lig = (LigandLight) ligands.get(u);

			LinkedList atomsL = lig.getAtoms();
			for (int i = 0; i < atomsL.size(); i++) {
				Atom atom = (Atom) atomsL.get(i);

				int x1 = getCx(atom.getFx());
				int x2 = x1 + 50;
				int y1 = getCy(atom.getFy());
				int y2 = y1 + 50;

				if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2)) {
					atomLigand = atom;
					return lig;
				}

			}

		}

		return null;
	}

	Atom getAtom(int x0, int y0) {

		int x1 = getCx(mainAtom.getFx());
		int x2 = x1 + 30;
		int y1 = getCy(mainAtom.getFy());
		int y2 = y1 + 10;
		if ((x0 >= x1) && (x0 <= x2) && (y0 >= y1) && (y0 <= y2))
			return mainAtom;
		return null;

	}

	public void rotation() {
		LinkedList all = getAll();
		LinkedList tmp = barycenter(all);
		rotation(all, (Double) tmp.get(0), (Double) tmp.get(1), -45);
		setScale();
		drawScale3();
		repaint();

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

	public LinkedList getAll() {

		LinkedList result = new LinkedList();

		result.add(mainAtom);

		for (int j = 0; j < ligands.size(); j++) {
			LigandLight lig = (LigandLight) ligands.get(j);

			LinkedList atomsE = lig.getAtoms();
			for (int i = 0; i < atomsE.size(); i++) {
				Atom ato = (Atom) atomsE.get(i);

				result.add(ato);
			}
		}

		for (int j = 0; j < water.size(); j++) {
			Atom w = (Atom) water.get(j);
			result.add(w);

		}

		for (int j = 0; j < resExplicit.size(); j++) {
			Explicit e = (Explicit) resExplicit.get(j);

			LinkedList atomsE = e.getAtoms();
			for (int i = 0; i < atomsE.size(); i++) {
				Atom atom = (Atom) atomsE.get(i);
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

}
