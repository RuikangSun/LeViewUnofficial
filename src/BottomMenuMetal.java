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
 This class defines the bottom menu for the metal frame;
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


public class BottomMenuMetal extends JPanel {

	Scrollbar contacts;
	JLabel l3;
	JLabel l4;
	JLabel currentC;
	JPanel boutons;
	JPanel boutons2;
	MyPanelMetal pan;

	public BottomMenuMetal(final MyFrameMetal f, final MyPanelMetal pan) {

		this.pan = pan;
		this.setBackground(Color.LIGHT_GRAY);
		setLayout(new GridLayout(3, 1, 10, 5));
		JLabel t2 = new JLabel("Close residue cutoff", JLabel.CENTER);
		add(t2);

		contacts = new Scrollbar(Scrollbar.HORIZONTAL, 0, 1, 0, 51);
		contacts.setValue(40);

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

		add(contacts);

		l3 = new JLabel();
		l3.setText("0");
		l3.setHorizontalTextPosition(JLabel.RIGHT);
		l4 = new JLabel();
		l4.setText("5.0");
		currentC = new JLabel("4.0", JLabel.CENTER);
		currentC.setForeground(Color.red);

		boutons = new JPanel(new BorderLayout());
		boutons.setBackground(Color.LIGHT_GRAY);

		boutons2 = new JPanel(new BorderLayout());
		boutons2.setBackground(Color.LIGHT_GRAY);

		boutons2.add(l3, BorderLayout.WEST);
		boutons2.add(currentC, BorderLayout.CENTER);
		boutons2.add(l4, BorderLayout.EAST);

		add(boutons2);

	}
}
