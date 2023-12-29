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
 This class defines the ligand frame.
 */

import java.awt.Color;
import java.io.*;
import java.util.*;

import javax.swing.JFrame;
import javax.swing.JPanel;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

public class MyFrame extends JFrame {
	public BufferedImage bi;
	public Ligand ligand;
	public String name;
	public MyPanel pan;
	public TopMenu panel;

	MainFrame mf;
	int index;

	public MyFrame(Ligand ligand) {

		this.ligand = ligand;
		mf=null;

		LinkedList residues = ligand.getResidues();
		for (int i = 0; i < residues.size(); i++) {
			Residue r = (Residue) residues.get(i);
			if (!r.getMapping()) {
				System.out.println("ATTENTION " + r.getName());

			}
		}

		Molecule mol = ligand.getMol();
		Chain c = ligand.getParent();
		
		name = c.getLongName();
		int pos = name.lastIndexOf("/");
		String name2;
		if (pos == -1)
			name2 = name;
		else
			name2 = name.substring(pos + 1, name.length());

		this.setTitle(name2);
		this.setSize(700, 700);
		this.setLocationRelativeTo(null);
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);


		pan = center();
		add(BorderLayout.CENTER, pan);

		buildLayout();

		this.repaint();

		
		this.setVisible(true);

	

	}

	public void buildLayout() {

		add(BorderLayout.NORTH, nord());
		add(BorderLayout.SOUTH, sud());

	}

	
	public MyPanel center() {
		MyPanel pan = new MyPanel(this, panel);
		getContentPane().add(pan.affiche(ligand));
		return pan;
	}

	public TopMenu nord() {
		TopMenu panel = new TopMenu(this, pan);
		return panel;
	}

	public BottomMenu sud() {
		BottomMenu panel = new BottomMenu(this, pan);
		return panel;
	}

	public void setMainFrame(MainFrame f) {
		mf = f;
		WindowListener listener = new WindowAdapter() {
				public void windowClosing(WindowEvent w) {
					mf.setCheckBox(index);
				}
			};
		
		this.addWindowListener(listener);
	}

	public void setIndex(int i) {
		index = i;
	}

}
