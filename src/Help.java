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
 This class defines the frame for the manual.
 */

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.tree.*;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import java.net.URL;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.net.*;
import java.io.*;
import java.io.File;

public class Help extends JFrame implements TreeSelectionListener {
	private JSplitPane split;
	JTree tree;
	JEditorPane editeur;
	JPanel pannel;
	JScrollPane pane;

	public Help() {

		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		DefaultMutableTreeNode root = new DefaultMutableTreeNode("Help");
		root.add(new DefaultMutableTreeNode("About the Method"));
		root.add(new DefaultMutableTreeNode("Options"));
		root.add(new DefaultMutableTreeNode("Water Mediated"));
		root.add(new DefaultMutableTreeNode("Changing cut-off"));
		root.add(new DefaultMutableTreeNode("Deleting interactions"));
		root.add(new DefaultMutableTreeNode("Moving Elements"));
		root.add(new DefaultMutableTreeNode("Rotation"));
		root.add(new DefaultMutableTreeNode("Export"));
		
		

		TreeModel model = new DefaultTreeModel(root);

		tree = new JTree(model);
		tree.addTreeSelectionListener(this);

		pane = new JScrollPane(tree);

		pannel = new JPanel();
		pannel.setBackground(Color.WHITE);

		try {

			URL monUrl = getClass().getResource("html/index.html");
			editeur = new JEditorPane(monUrl);

			editeur.setEditable(false);
			editeur.addHyperlinkListener(new HyperlinkListener() {
				public void hyperlinkUpdate(HyperlinkEvent e) {
					if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
						URL url = e.getURL();
						if (url == null)
							return;
						try {
							editeur.setPage(e.getURL());
						} catch (Exception ex) {
							ex.printStackTrace();
						}
					}
				}
			});

			pannel.add(editeur);
		} catch (Exception e1) {
			e1.printStackTrace();
		}

		split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, pane,
				new JScrollPane(pannel));

		this.getContentPane().add(split, BorderLayout.CENTER);

		setSize(650, 650);
		setVisible(true);
	}

	public void valueChanged(TreeSelectionEvent event) {
		String s = tree.getLastSelectedPathComponent().toString();
		if (s.equals("Export")) {

			try {
				URL monUrl = getClass().getResource("html/save.html");
				editeur.setPage(monUrl);
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			pannel.repaint();

		}

		if (s.equals("Help")) {

			try {
				URL monUrl = getClass().getResource("html/index.html");
				editeur.setPage(monUrl);
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			pannel.repaint();

		}

		if (s.equals("Options")) {

			try {
				URL monUrl = getClass().getResource("html/options.html");
				editeur.setPage(monUrl);
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			pannel.repaint();

		}

		if (s.equals("Deleting interactions")) {

			try {
				URL monUrl = getClass().getResource("html/delete.html");
				editeur.setPage(monUrl);
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			pannel.repaint();

		}

		if (s.equals("About the Method")) {

			try {
				URL monUrl = getClass().getResource("html/about.html");
				editeur.setPage(monUrl);
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			pannel.repaint();

		}

		if (s.equals("Changing cut-off")) {

			try {
				URL monUrl = getClass().getResource("html/cutoff.html");
				editeur.setPage(monUrl);
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			pannel.repaint();

		}

		if (s.equals("Moving Elements")) {

			try {
				URL monUrl = getClass().getResource("html/moving.html");
				editeur.setPage(monUrl);
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			pannel.repaint();

		}

		if (s.equals("Water Mediated")) {

			try {
				URL monUrl = getClass().getResource("html/waterMediated.html");
				editeur.setPage(monUrl);
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			pannel.repaint();

		}

		if (s.equals("Rotation")) {

			try {
				URL monUrl = getClass().getResource("html/rotation.html");
				editeur.setPage(monUrl);
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			pannel.repaint();

		}

	}

}
