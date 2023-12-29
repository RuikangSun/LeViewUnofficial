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
 This class defines the Bottom menu with the cutoff sliders.
 */

import java.awt.event.*;
import javax.swing.*;
import java.awt.*;
import java.awt.Color;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Component;
import javax.swing.JPanel;

import java.io.*;
import java.util.*;

public class BottomMenu extends JPanel {

	Scrollbar Hbonds;
	Scrollbar contacts;
	JLabel l1;
	JLabel l2;
	JLabel l3;
	JLabel l4;
	JLabel currentHB;
	JLabel currentC;
	JPanel boutons;
	JPanel boutons2;
	MyPanel pan;

	public BottomMenu(final MyFrame f, final MyPanel pan) {

		this.pan = pan;
		this.setBackground(Color.LIGHT_GRAY);
		setLayout(new GridLayout(3, 2, 10, 5));

		JLabel t1 = new JLabel("Hbonds cutoff", JLabel.CENTER);
		add(t1);
		JLabel t2 = new JLabel("Close residue cutoff", JLabel.CENTER);
		add(t2);

		Hbonds = new Scrollbar(Scrollbar.HORIZONTAL, 0, 1, 0, 51);
		Hbonds.setValue(33);

		/*
		 * if the mode MOVE ELEMENTS is active the sliders is blocked on the
		 * maximal value
		 */
		AdjustmentListener hListener = new AdjustmentListener() {
			public void adjustmentValueChanged(AdjustmentEvent e) {
				if (pan.getLockedMode() == 1) {
					int v = (int) (pan.getValueHB() * 10);
					if (Hbonds.getValue() > v) {
						Hbonds.setValue(v);
					} else {
						double temp = Hbonds.getValue() / 10.0;
						currentHB.setText(String.valueOf(temp));
						pan.setCutoffHB(temp);
						pan.repaint();
					}
				} else {
					double temp = Hbonds.getValue() / 10.0;
					currentHB.setText(String.valueOf(temp));
					pan.setCutoffHB(temp);
					pan.repaint();
				}

			}
		};

		Hbonds.addAdjustmentListener(hListener);

		contacts = new Scrollbar(Scrollbar.HORIZONTAL, 0, 1, 0, 51);
		contacts.setValue(40);

		/*
		 * if the mode MOVE ELEMENTS is active the sliders is blocked on the
		 * maximal value
		 */
		AdjustmentListener hListener2 = new AdjustmentListener() {
			public void adjustmentValueChanged(AdjustmentEvent e) {
				if (pan.getLockedMode() == 1) {
					int v = (int) (pan.getValueC() * 10);
					if (contacts.getValue() > v) {
						contacts.setValue(v);
					} else {
						double temp = contacts.getValue() / 10.0;
						currentC.setText(String.valueOf(temp));
						pan.setCutoffC(temp);
						pan.repaint();
					}
				} else {
					double temp = contacts.getValue() / 10.0;
					currentC.setText(String.valueOf(temp));
					pan.setCutoffC(temp);
					pan.repaint();
				}

			}
		};

		contacts.addAdjustmentListener(hListener2);

		add(Hbonds);
		add(contacts);

		l1 = new JLabel();
		l1.setText("0");
		l2 = new JLabel();
		l2.setText("5.0");
		currentHB = new JLabel("3.3", JLabel.CENTER);
		currentHB.setForeground(Color.red);
		l3 = new JLabel();
		l3.setText("0");
		l3.setHorizontalTextPosition(JLabel.RIGHT);
		l4 = new JLabel();
		l4.setText("5.0");
		currentC = new JLabel("4.0", JLabel.CENTER);
		currentC.setForeground(Color.red);

		boutons = new JPanel(new BorderLayout());
		boutons.setBackground(Color.LIGHT_GRAY);
		boutons.add(l1, BorderLayout.WEST);
		boutons.add(currentHB, BorderLayout.CENTER);
		boutons.add(l2, BorderLayout.EAST);

		add(boutons);

		boutons2 = new JPanel(new BorderLayout());
		boutons2.setBackground(Color.LIGHT_GRAY);

		boutons2.add(l3, BorderLayout.WEST);
		boutons2.add(currentC, BorderLayout.CENTER);
		boutons2.add(l4, BorderLayout.EAST);

		add(boutons2);
	}
}
