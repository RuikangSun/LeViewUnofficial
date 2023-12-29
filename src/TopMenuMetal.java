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
 This class defines the top menu for the metal frame.
 */

import java.awt.event.*;
import javax.swing.*;
import java.awt.*;

import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.JPanel;

import java.io.*;
import java.util.*;

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

import java.awt.Dimension;
import java.io.File;
import java.util.Properties;
import org.freehep.graphics2d.VectorGraphics;
import org.freehep.util.export.ExportDialog;

import java.text.*;
public class TopMenuMetal extends JPanel {

	MyPanelMetal pan;

	String name;

	int mode;
	Color curMetal;
	Color curEx;
	Color curRes;
	Color curName;
	Color curName2;
	Color curLig;
	Color brown;
	Color curWaters;

	public TopMenuMetal(final MyFrameMetal f, final MyPanelMetal pan) {

		name = f.name;
		this.pan = pan;

		brown = new Color(139, 71, 38);

		BorderLayout gestionnaire = new BorderLayout();
		this.setLayout(gestionnaire);

		curMetal = Color.GREEN;
		curEx = Color.BLACK;
		curLig = Color.GRAY;
		curRes = Color.WHITE;
		curName = Color.GRAY;
		curName2 = Color.BLACK;
		curWaters = brown;

		JMenuBar menuBar = new JMenuBar();

		JMenu save = new JMenu("Export");
		JMenu colour = new JMenu("Options");
		JMenu color = new JMenu("Residue Colour Scheme");
		JMenu help = new JMenu("Help");

		JMenu labels = new JMenu("Atom Labels");
		JMenu atomColour = new JMenu("Atom Colour Scheme");

		JMenuItem about = new JMenuItem("About");
		JMenuItem manual = new JMenuItem("Manual");

		JMenu change = new JMenu("Element Colour");
		JMenu style = new JMenu("HBond style");

		JMenuItem metal = new JMenuItem("Metal");
		JMenuItem resExplicit = new JMenuItem("Explicit residues");
		JMenuItem ligands = new JMenuItem("Ligands");
		JMenuItem waters = new JMenuItem("Waters");
		JMenuItem res = new JMenuItem("Close residues");
		JMenuItem name = new JMenuItem("Ligand names");
		JMenuItem name2 = new JMenuItem("Explicit residue names");
		change.add(metal);
		change.add(resExplicit);
		change.add(ligands);
		change.add(waters);
		change.add(res);
		change.add(name);
		change.add(name2);

		ActionListener labelActionListener = new ActionListener() {
			public void actionPerformed(ActionEvent actionEvent) {
				AbstractButton aButton = (AbstractButton) actionEvent
						.getSource();
				String s = aButton.getText();
				if (s.equals("Display")) {
					pan.setLabelMode(0);
					pan.repaint();
				}
				if (s.equals("Hide")) {
					pan.setLabelMode(1);
					pan.repaint();
				}
			}
		};

		JRadioButtonMenuItem item2;
		ButtonGroup labelGroup = new ButtonGroup();
		item2 = new JRadioButtonMenuItem("Hide", true);
		labels.add(item2);
		labelGroup.add(item2);
		item2.addActionListener(labelActionListener);
		item2 = new JRadioButtonMenuItem("Display");
		labels.add(item2);
		labelGroup.add(item2);
		item2.addActionListener(labelActionListener);

		ActionListener atomActionListener = new ActionListener() {
			public void actionPerformed(ActionEvent actionEvent) {
				AbstractButton aButton = (AbstractButton) actionEvent
						.getSource();
				String s = aButton.getText();
				if (s.equals("Standard")) {
					pan.setAtomColorMode(0);
					pan.repaint();
				}
				if (s.equals("Plain Colour")) {
					pan.setAtomColorMode(1);
					pan.repaint();
				}
			}
		};

		JRadioButtonMenuItem item3;
		ButtonGroup atomGroup = new ButtonGroup();
		item3 = new JRadioButtonMenuItem("Standard", true);
		atomColour.add(item3);
		atomGroup.add(item3);
		item3.addActionListener(atomActionListener);
		item3 = new JRadioButtonMenuItem("Plain Colour");
		atomColour.add(item3);
		atomGroup.add(item3);
		item3.addActionListener(atomActionListener);

		metal.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Color newColor = JColorChooser.showDialog(pan,
						"Choose Metal Colour", curMetal);
				if (newColor != null) {
					pan.setMetalColor(newColor);
					pan.repaint();
					curMetal = newColor;
				}
			}
		});

		resExplicit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Color newColor = JColorChooser.showDialog(pan,
						"Choose Explicit residue Colour", curEx);
				if (newColor != null) {
					pan.setResExplicitColor(newColor);
					pan.repaint();
					curEx = newColor;
				}
			}
		});

		ligands.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Color newColor = JColorChooser.showDialog(pan,
						"Choose Ligand Colour", curLig);
				if (newColor != null) {
					pan.setLigandsColor(newColor);
					pan.repaint();
					curLig = newColor;
				}
			}
		});

		waters.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Color newColor = JColorChooser.showDialog(pan,
						"Choose water Colour", curWaters);
				if (newColor != null) {
					pan.setWatersColor(newColor);
					pan.repaint();
					curWaters = newColor;
				}
			}
		});

		res.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Color newColor = JColorChooser.showDialog(pan,
						"Choose close residue Colour", curRes);
				if (newColor != null) {
					pan.setColorMode(5);
					pan.setResColor(newColor);
					pan.repaint();
					curRes = newColor;
				}
			}
		});

		name.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Color newColor = JColorChooser.showDialog(pan,
						"Choose ligand name Colour", curName);
				if (newColor != null) {
					pan.setNameColor(newColor);
					pan.repaint();
					curName = newColor;
				}
			}
		});

		name2.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Color newColor = JColorChooser.showDialog(pan,
						"Choose explicit name Colour", curName2);
				if (newColor != null) {
					pan.setName2Color(newColor);
					pan.repaint();
					curName2 = newColor;
				}
			}
		});

		about.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveImagePNG();

			}
		});

		manual.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Help h = new Help();

			}
		});

		help.add(manual);

		JMenuItem png = new JMenuItem(".png");
		JMenuItem jpeg = new JMenuItem(".jpg");
		JMenuItem gif = new JMenuItem(".gif");
		JMenuItem pdf = new JMenuItem(".pdf");
		JMenuItem svg = new JMenuItem(".svg");
		JMenuItem eps = new JMenuItem(".eps");
		JMenuItem flat = new JMenuItem("flat");

		JMenuItem hydro = new JMenuItem("Hydro.");
		JMenuItem prop = new JMenuItem("Prop.");
		JMenuItem charge = new JMenuItem("Charge");
		JMenuItem struct = new JMenuItem("Struct.");

		save.add(png);
		save.add(jpeg);
		save.add(gif);
		save.add(pdf);
		save.add(svg);
		save.add(eps);
		save.add(flat);

		color.add(hydro);
		color.add(prop);
		color.add(charge);
		color.add(struct);

		colour.add(labels);
		colour.add(atomColour);
		colour.add(color);
		colour.add(change);

		menuBar.add(save);
		menuBar.add(colour);
		menuBar.add(help);

		add(menuBar, BorderLayout.WEST);

		final JToggleButton modif = new JToggleButton("Move elements");

		add(modif, BorderLayout.CENTER);

		ImageIcon icon = new ImageIcon(this.getClass().getResource(
				"rotation.png"));
		JButton rotation = new JButton(icon);

		add(rotation, BorderLayout.EAST);

		png.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveImagePNG();

			}
		});

		rotation.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				pan.rotation();

			}
		});

		pdf.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveImagePDF();

			}
		});

		svg.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveImageSVG();

			}
		});

		eps.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveImageEPS();

			}
		});


		flat.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveImageFlat();

			}
		});

		modif.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (modif.getLabel().equals("Move elements")) {
					modif.setLabel("Cancel");
					pan.setLockedMode(1);
				} else {
					int choix = JOptionPane.showConfirmDialog(f,
							"Are you sure? You may lose your modifications.",
							"Confirm", JOptionPane.OK_CANCEL_OPTION,
							JOptionPane.QUESTION_MESSAGE);

					if (choix == 0) {

						modif.setLabel("Move elements");
						pan.setLockedMode(0);
					}

				}
			}
		});

		jpeg.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveImageJPEG();

			}
		});

		gif.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveImageGIF();

			}
		});

		hydro.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				pan.setColorMode(1);
				pan.repaint();

			}
		});

		prop.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				pan.setColorMode(2);
				pan.repaint();

			}
		});

		charge.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				pan.setColorMode(3);
				pan.repaint();

			}
		});

		struct.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				pan.setColorMode(4);
				pan.repaint();

			}
		});

	}

	public boolean saveImagePNG() {
		File fileOut = new File(name + ".png");

		JFileChooser fOut = new JFileChooser();

		fOut.setApproveButtonText("Select");
		fOut.setSelectedFile(fileOut);
		int button = fOut.showSaveDialog(new JFrame());

		while (button == JFileChooser.APPROVE_OPTION) {
			fileOut = fOut.getSelectedFile();

			if (fileOut.getName().indexOf(".") == -1) {
				fileOut.renameTo(new File(fileOut.getAbsolutePath() + ".png"));
			}

			if (!fileOut.exists()) {
				try {
					fileOut.createNewFile();
				} catch (IOException e) {
					JOptionPane.showMessageDialog(null, "Error creating file",
							"Warning", JOptionPane.ERROR_MESSAGE);
					return false;
				}

				save(pan, fileOut);
				return true;
			} else {

				int confirm = JOptionPane.showConfirmDialog(null,
						"This file already exist, do you want to erase it ?");
				if (confirm == 0) {

					save(pan, fileOut);
					return true;
				} else {
					JOptionPane.showMessageDialog(null,
							"Please select another file", "Warning",
							JOptionPane.ERROR_MESSAGE);
				}
			}
			fOut.setSelectedFile(fileOut);
			button = fOut.showOpenDialog(new JFrame());
		}
		return false;
	}

	public boolean saveImageJPEG() {
		File fileOut = new File(name + ".jpg");
		JFileChooser fOut = new JFileChooser();

		fOut.setApproveButtonText("Select");
		fOut.setSelectedFile(fileOut);
		int button = fOut.showSaveDialog(new JFrame());

		while (button == JFileChooser.APPROVE_OPTION) {
			fileOut = fOut.getSelectedFile();

			if (fileOut.getName().indexOf(".") == -1) {
				fileOut.renameTo(new File(fileOut.getAbsolutePath() + ".jpg"));
			}

			if (!fileOut.exists()) {
				try {
					fileOut.createNewFile();
				} catch (IOException e) {
					JOptionPane.showMessageDialog(null, "Error creating file",
							"Warning", JOptionPane.ERROR_MESSAGE);
					return false;
				}

				save2(pan, fileOut);
				return true;
			} else {

				int confirm = JOptionPane.showConfirmDialog(null,
						"This file already exist, do you want to erase it ?");
				if (confirm == 0) {

					save2(pan, fileOut);
					return true;
				} else {
					JOptionPane.showMessageDialog(null,
							"Please select another file", "Warning",
							JOptionPane.ERROR_MESSAGE);
				}
			}
			fOut.setSelectedFile(fileOut);
			button = fOut.showOpenDialog(new JFrame());
		}
		return false;
	}

	public boolean saveImageGIF() {
		File fileOut = new File(name + ".gif");
		JFileChooser fOut = new JFileChooser();

		fOut.setApproveButtonText("Select");
		fOut.setSelectedFile(fileOut);
		int button = fOut.showSaveDialog(new JFrame());

		while (button == JFileChooser.APPROVE_OPTION) {
			fileOut = fOut.getSelectedFile();

			if (fileOut.getName().indexOf(".") == -1) {
				fileOut.renameTo(new File(fileOut.getAbsolutePath() + ".gif"));
			}

			if (!fileOut.exists()) {
				try {
					fileOut.createNewFile();
				} catch (IOException e) {
					JOptionPane.showMessageDialog(null, "Error creating file",
							"Warning", JOptionPane.ERROR_MESSAGE);
					return false;
				}

				save3(pan, fileOut);
				return true;
			} else {
				int confirm = JOptionPane.showConfirmDialog(null,
						"This file already exist, do you want to erase it ?");
				if (confirm == 0) {

					save3(pan, fileOut);
					return true;
				} else {
					JOptionPane.showMessageDialog(null,
							"Please select another file", "Warning",
							JOptionPane.ERROR_MESSAGE);
				}
			}
			fOut.setSelectedFile(fileOut);
			button = fOut.showOpenDialog(new JFrame());
		}
		return false;
	}

	public boolean saveImagePDF() {
		File fileOut = new File(name + ".pdf");

		JFileChooser fOut = new JFileChooser();

		fOut.setApproveButtonText("Select");
		fOut.setSelectedFile(fileOut);
		int button = fOut.showSaveDialog(new JFrame());

		while (button == JFileChooser.APPROVE_OPTION) {
			fileOut = fOut.getSelectedFile();

			if (fileOut.getName().indexOf(".") == -1) {
				fileOut.renameTo(new File(fileOut.getAbsolutePath() + ".pdf"));
			}

			if (!fileOut.exists()) {
				try {
					fileOut.createNewFile();
				} catch (IOException e) {
					JOptionPane.showMessageDialog(null, "Error creating file",
							"Warning", JOptionPane.ERROR_MESSAGE);
					return false;
				}

				save4(pan, fileOut);
				return true;
			} else {

				int confirm = JOptionPane.showConfirmDialog(null,
						"This file already exist, do you want to erase it ?");
				if (confirm == 0) {

					save4(pan, fileOut);
					return true;
				} else {
					JOptionPane.showMessageDialog(null,
							"Please select another file", "Warning",
							JOptionPane.ERROR_MESSAGE);
				}
			}
			fOut.setSelectedFile(fileOut);
			button = fOut.showOpenDialog(new JFrame());
		}
		return false;
	}

	public boolean saveImageSVG() {
		File fileOut = new File(name + ".svg");
		JFileChooser fOut = new JFileChooser();

		fOut.setApproveButtonText("Select");
		fOut.setSelectedFile(fileOut);
		int button = fOut.showSaveDialog(new JFrame());

		while (button == JFileChooser.APPROVE_OPTION) {
			fileOut = fOut.getSelectedFile();

			if (fileOut.getName().indexOf(".") == -1) {
				fileOut.renameTo(new File(fileOut.getAbsolutePath() + ".svg"));
			}

			if (!fileOut.exists()) {
				try {
					fileOut.createNewFile();
				} catch (IOException e) {
					JOptionPane.showMessageDialog(null, "Error creating file",
							"Warning", JOptionPane.ERROR_MESSAGE);
					return false;
				}

				save5(pan, fileOut);
				return true;
			} else {
				int confirm = JOptionPane.showConfirmDialog(null,
						"This file already exist, do you want to erase it ?");
				if (confirm == 0) {

					save5(pan, fileOut);
					return true;
				} else {
					JOptionPane.showMessageDialog(null,
							"Please select another file", "Warning",
							JOptionPane.ERROR_MESSAGE);
				}
			}
			fOut.setSelectedFile(fileOut);
			button = fOut.showOpenDialog(new JFrame());
		}
		return false;
	}

	public boolean saveImageEPS() {
		File fileOut = new File(name + ".eps");
		JFileChooser fOut = new JFileChooser();

		fOut.setApproveButtonText("Select");
		fOut.setSelectedFile(fileOut);
		int button = fOut.showSaveDialog(new JFrame());

		while (button == JFileChooser.APPROVE_OPTION) {
			fileOut = fOut.getSelectedFile();

			if (fileOut.getName().indexOf(".") == -1) {
				fileOut.renameTo(new File(fileOut.getAbsolutePath() + ".eps"));
			}

			if (!fileOut.exists()) {
				try {
					fileOut.createNewFile();
				} catch (IOException e) {
					JOptionPane.showMessageDialog(null, "Error creating file",
							"Warning", JOptionPane.ERROR_MESSAGE);
					return false;
				}

				save6(pan, fileOut);
				return true;
			} else {
				int confirm = JOptionPane.showConfirmDialog(null,
						"This file already exist, do you want to erase it ?");
				if (confirm == 0) {

					save6(pan, fileOut);
					return true;
				} else {
					JOptionPane.showMessageDialog(null,
							"Please select another file", "Warning",
							JOptionPane.ERROR_MESSAGE);
				}
			}
			fOut.setSelectedFile(fileOut);
			button = fOut.showOpenDialog(new JFrame());
		}
		return false;
	}

	public void save(JPanel jp, File fic) {

		BufferedImage bufferedImage;
		bufferedImage = new BufferedImage(jp.getWidth(), jp.getHeight(),
				BufferedImage.TYPE_INT_RGB);

		Graphics g = bufferedImage.createGraphics();
		jp.paint(g);

		try {
			ImageIO.write(bufferedImage, "png", fic);
		}

		catch (Exception e) {
			System.out.println("erreur enregistrement image...");
			e.printStackTrace();
		}

	}



	public boolean saveImageFlat() {
		File fileOut = new File(name + ".txt");

		JFileChooser fOut = new JFileChooser();

		fOut.setApproveButtonText("Select");
		fOut.setSelectedFile(fileOut);
		int button = fOut.showSaveDialog(new JFrame());

		while (button == JFileChooser.APPROVE_OPTION) {
			fileOut = fOut.getSelectedFile();

			if (fileOut.getName().indexOf(".") == -1) {
				fileOut.renameTo(new File(fileOut.getAbsolutePath() + ".txt"));
			}

			if (!fileOut.exists()) {
				try {
					fileOut.createNewFile();
				} catch (IOException e) {
					JOptionPane.showMessageDialog(null, "Error creating file",
							"Warning", JOptionPane.ERROR_MESSAGE);
					return false;
				}

				saveFlat(pan, fileOut.toString());
				return true;
			} else {

				int confirm = JOptionPane.showConfirmDialog(null,
						"This file already exist, do you want to erase it ?");
				if (confirm == 0) {

					saveFlat(pan, fileOut.toString());
					return true;
				} else {
					JOptionPane.showMessageDialog(null,
							"Please select another file", "Warning",
							JOptionPane.ERROR_MESSAGE);
				}
			}
			fOut.setSelectedFile(fileOut);
			button = fOut.showOpenDialog(new JFrame());
		}
		return false;
	}





	public void saveFlat(MyPanelMetal pan,String fic) {

		try {
		PrintWriter out=new PrintWriter(new FileWriter(fic));

		

		LinkedList cl=pan.explicit;
		if(cl.size()!=0){
		out.println("*****************Ligand  atoms*****************");
		out.println("metal_atom\tligand_atom\tresidue\tdistance");
		for(int i=0;i<cl.size();i++){
			Bond b=(Bond)cl.get(i);

			Atom don=b.getFirst();
			Atom acc=b.getSecond();
			
			Residue res_acc=acc.getParent();
			Chain c_don=res_acc.getParent();
			DecimalFormat df = new DecimalFormat("########.0000");
			String distance = df.format((don.distance(acc)));
			Chain c_acc=res_acc.getParent();

			out.println(don.getName()+"_"+don.getId()+"("+c_don.getName()+")"+
			"\t"+acc.getName()+"_"+acc.getId()+"\t"+res_acc.getName()+"_"+res_acc.getId()+"("+c_don.getName()+")"+"\t"+distance);

		}
		}


		 cl=pan.contactList;

		LinkedList tmp=new LinkedList();
		for(int i=0;i<cl.size();i++){
			Contact c=(Contact)cl.get(i);
			double dd = c.getDistance() / 2;
			
			if (dd <= pan.cutoffC) {
				tmp.add(c);

			}
		}


		if(tmp.size()!=0){
		out.println("*****************Nearby  residues*****************");
		out.println("metal_atom\tresidue\tnearby_atom\tresidue\tdistance");

		for(int i=0;i<tmp.size();i++){
			Contact c=(Contact)tmp.get(i);
			Atom don=c.getRef();
			Atom acc=c.getTrait();
			Residue res_don=don.getParent();
			Residue res_acc=acc.getParent();
			Chain c_don=res_acc.getParent();
			DecimalFormat df = new DecimalFormat("########.0000");
			String distance = df.format((c.getDistance()/2));
			Chain c_acc=res_acc.getParent();

			out.println(don.getName()+"_"+don.getId()+"\t"+res_don.getName()+"_"+res_don.getId()+"("+c_don.getName()+")"+
			"\t"+acc.getName()+"_"+acc.getId()+"\t"+res_acc.getName()+"_"+res_acc.getId()+"("+c_don.getName()+")"+"\t"+distance);
			
		}
		}



		 cl=pan.explicit2;
		if(cl.size()!=0){
		out.println("*****************H2O  atoms*****************");
		out.println("metal_atom\tH2O_atom\tresidue\tdistance");
		for(int i=0;i<cl.size();i++){
			Bond b=(Bond)cl.get(i);

			Atom don=b.getFirst();
			Atom acc=b.getSecond();
			
			//Residue res_acc=acc.getParent();
			//Chain c_don=res_acc.getParent();
			DecimalFormat df = new DecimalFormat("########.0000");
			String distance = df.format((don.distance(acc)));
			//Chain c_acc=res_acc.getParent();

			out.println(don.getName()+"_"+don.getId()+
			"\t"+acc.getName()+"_"+acc.getId()+"\t"+distance);

		}
		}


		out.close();
		}
		catch (Exception e) {
			System.out.println("erreur flat format...");
			e.printStackTrace();
		}
		

	}


	public void save2(JPanel jp, File fic) {

		BufferedImage bufferedImage;
		bufferedImage = new BufferedImage(jp.getWidth(), jp.getHeight(),
				BufferedImage.TYPE_INT_RGB);

		Graphics g = bufferedImage.createGraphics();
		jp.paint(g);

		try {
			ImageIO.write(bufferedImage, "JPEG", fic);
		}

		catch (Exception e) {
			System.out.println("erreur enregistrement image...");
			e.printStackTrace();
		}

	}

	public void save3(JPanel jp, File fic) {

		BufferedImage bufferedImage;
		bufferedImage = new BufferedImage(jp.getWidth(), jp.getHeight(),
				BufferedImage.TYPE_INT_RGB);

		Graphics g = bufferedImage.createGraphics();
		jp.paint(g);

		try {
			ImageIO.write(bufferedImage, "GIF", fic);
		}

		catch (Exception e) {
			System.out.println("erreur enregistrement image...");
			e.printStackTrace();
		}

	}

	public void save4(JPanel jp, File fic) {

		try {
			Properties p = new Properties();

			VectorGraphics g = new org.freehep.graphicsio.pdf.PDFGraphics2D(
					fic, new Dimension(jp.size().width, jp.size().height));

			g.startExport();
			jp.print(g);
			g.endExport();
		} catch (Exception ex) {
			System.out.println("erreur enregistrement image...");
			ex.printStackTrace();
		}

	}

	public void save5(JPanel jp, File fic) {

		try {
			Properties p = new Properties();
			p.setProperty("PageSize", "A5");
			VectorGraphics g = new org.freehep.graphicsio.svg.SVGGraphics2D(
					fic, new Dimension(jp.size().width, jp.size().height));
			g.setProperties(p);
			g.startExport();
			jp.print(g);
			g.endExport();
		} catch (Exception ex) {
			System.out.println("erreur enregistrement image...");
			ex.printStackTrace();
		}

	}

	public void save6(JPanel jp, File fic) {

		try {
			Properties p = new Properties();
			p.setProperty("PageSize", "A5");
			VectorGraphics g = new org.freehep.graphicsio.ps.PSGraphics2D(fic,
					new Dimension(jp.size().width, jp.size().height));

			g.setProperties(p);
			g.startExport();
			jp.print(g);
			g.endExport();
		} catch (Exception ex) {
			System.out.println("erreur enregistrement image...");
			ex.printStackTrace();
		}

	}

}
