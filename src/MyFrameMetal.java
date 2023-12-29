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
 This class defines the metal frame.
 */

import java.awt.Color;

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

public class MyFrameMetal extends JFrame {
	public BufferedImage bi;
	public Metal metal;
	public String name;
	public MyPanelMetal pan;
	public TopMenuMetal panel;

	MainFrame mf;
	int index;

	public MyFrameMetal(Metal metal) {

		this.metal = metal;
		Molecule mol = metal.getMol();
		Chain c = metal.getParent();
		name = c.getLongName();
		int pos = name.lastIndexOf("/");
		String name2;
		if (pos == -1)
			name2 = name;
		else
			name2 = name.substring(pos + 1, name.length());

		this.setTitle(name2);
		this.setSize(600, 600);
		this.setLocationRelativeTo(null);
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		WindowListener listener = new WindowAdapter() {
			public void windowClosing(WindowEvent w) {
				mf.setCheckBox2(index);
			}
		};
		this.addWindowListener(listener);

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

	public MyPanelMetal center() {
		MyPanelMetal pan = new MyPanelMetal(this, panel);
		getContentPane().add(pan.affiche(metal));
		return pan;
	}

	public TopMenuMetal nord() {
		TopMenuMetal panel = new TopMenuMetal(this, pan);
		return panel;
	}

	public BottomMenuMetal sud() {
		BottomMenuMetal panel = new BottomMenuMetal(this, pan);
		return panel;
	}

	public void setMainFrame(MainFrame f) {
		mf = f;
	}

	public void setIndex(int i) {
		index = i;
	}

}
