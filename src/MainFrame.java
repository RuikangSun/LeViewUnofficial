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
 This class defines the main frame which displays all the ligands and ions.
 */

import javax.swing.*;
import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.BorderFactory;
import javax.swing.JProgressBar;
import javax.swing.border.Border;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.beans.*;
import java.util.Random;
 

public class MainFrame extends JFrame {

	 JCheckBox[] tab;
	 JCheckBox[] tab2;
	 MyFrame[] tabFl;
	 MyFrameMetal[] tabFm;
	JScrollPane scrollPane;
	JPanel pane;
	
	MainFrame frame;
	

	
	public MainFrame(String fic) throws IOException {

		frame = this;
		LinkedList ligands = new LinkedList();
		LinkedList metals = new LinkedList();
		
		pane = new JPanel();

		

		File f = new File(fic);
		Molecule m = f.getMolecule();
		
		Traitement t = new Traitement(m);
		

		LinkedList chains = m.getChains();
		for (int i = 0; i < chains.size(); i++) {
			Chain c = (Chain) chains.get(i);

			if (c.getType().equals("ligand")) {

				

				ligands.add(c);

				
			}

			if (c.getType().equals("ion")) {

				

				metals.add(c);

				
			}
		}

		

		int ligne = ligands.size() + metals.size() + 3;

		pane.setLayout(new GridLayout(ligne, 1));
		
		String nameMol="";
		if(m.getName().indexOf("/")!=-1){
			String ttt=m.getName();
			String [] tt=ttt.split("/");
			nameMol=tt[tt.length-1];
		
		}else{

			if(m.getName().indexOf("\\")!=-1){
			String ttt=m.getName();
			String [] tt=ttt.split("\\\\");
			nameMol=tt[tt.length-1];
		
			}else{
			nameMol=m.getName();}
		}
		JLabel title = new JLabel("PDB ID: " + nameMol, JLabel.CENTER);
		pane.add(title);

		JLabel l1;

		if (ligands.size() == 1) {
			l1 = new JLabel(ligands.size() + " ligand identified");
		} else {
			if (ligands.size() == 0) {
				l1 = new JLabel("no ligand");
			} else {
				l1 = new JLabel(ligands.size() + " ligands identified");
			}

		}
		pane.add(l1);

		tabFl = new MyFrame[ligands.size()];
		tab = new JCheckBox[ligands.size()];
		for (int i = 0; i < ligands.size(); i++) {

			
			final Chain c = (Chain) ligands.get(i);
			
			String name = c.getLongName2();
			tab[i] = new JCheckBox(name);
			
			final int ind = i;
			final String sego=m.getName();
			tab[i].addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					
					if (tab[ind].isSelected()) {

						final Ligand ligand = new Ligand(c, 4.0, 3.3);

						setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
						


						

						MyFrame fl = new MyFrame(ligand);
						
						setCursor(null);
						
						
						
						

						fl.setMainFrame(frame);
						fl.setIndex(ind);
						tabFl[ind] = fl;

						




					} else {
						if (tabFl[ind] != null) {
							tabFl[ind].dispose();
						}
					}

				}
			});
			pane.add(tab[i]);
		}

		
		JLabel l2;

		if (metals.size() == 1) {
			l2 = new JLabel(metals.size() + " ion identified");
		} else {
			if (metals.size() == 0) {
				l2 = new JLabel("no ion");
			} else {
				l2 = new JLabel(metals.size() + " ions identified");
			}

		}
		pane.add(l2);

		tabFm = new MyFrameMetal[metals.size()];
		tab2 = new JCheckBox[metals.size()];
		for (int i = 0; i < metals.size(); i++) {

			
			final Chain c = (Chain) metals.get(i);
			
			String name = c.getLongName2();
			tab2[i] = new JCheckBox(name);
			
			final int ind = i;
			tab2[i].addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					if (tab2[ind].isSelected()) {
						final Metal metal = new Metal(c, 4.0);
						setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
						MyFrameMetal fm = new MyFrameMetal(metal);
						setCursor(null);
						fm.setMainFrame(frame);
						fm.setIndex(ind);
						tabFm[ind] = fm;
					} else {
						if (tabFm[ind] != null) {
							tabFm[ind].dispose();
						}
					}

				}
			});
			pane.add(tab2[i]);
		}
		
		scrollPane = new JScrollPane(pane);
		this.getContentPane().add(scrollPane);
		
		this.setTitle("LeView");
		
		this.setSize(400, 400);
		//this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		this.setVisible(true);

	}

	public void setCheckBox(int i) {

		tab[i].setSelected(false);
		tabFl[i] = null;
	}

	public void setCheckBox2(int i) {

		tab2[i].setSelected(false);
		tabFm[i] = null;
	}


	

}
