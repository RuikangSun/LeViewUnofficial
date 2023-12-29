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

/* This is the main frame of LeView */

import javax.swing.*;
import java.awt.event.*;
import java.awt.*;
import java.io.*;
import java.io.File;
import javax.swing.filechooser.*;
import java.net.URL;
import java.util.*;
import javax.imageio.ImageIO;
import java.awt.Graphics;
import java.awt.image.BufferedImage;

public class LeView extends JFrame implements ActionListener
{
   JTextField fichier; 
	JTextField iden; 
   JButton search;   
   JButton exit;
   JButton open;   
   JButton download;   
   JFileChooser chooser;
   String fileName; 

	JTextArea log; 
	JScrollPane sbrText;

	String com;

	JRadioButton rb;
	JTextField res1;
	JTextField res2;
	 JTextField chain;
	
 
   public LeView ()
   {

	
      super("LeView: Ligand Environment Viewer");
      this.setSize(570,460);
      this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	fileName="";
	com="";
 
      BorderLayout agencement = new BorderLayout();
      this.setLayout(agencement);
 
      JPanel haut = new JPanel();
	JPanel milieu = new JPanel();
      JPanel bas = new JPanel();

	JPanel line1=new JPanel();
	JPanel line2=new JPanel();
	JPanel line3=new JPanel();

	
 
      /* on utilise un FlowLayout pour mettre les trois composants
       * les uns a la suite des autres (JLabel, JTextField, et JButton)
       */
      FlowLayout layout = new FlowLayout();
	GridLayout gl=new GridLayout(3,1);
      	haut.setLayout(layout);
	 milieu.setLayout(gl);
	bas.setLayout(layout);

	milieu.add(line1);
	milieu.add(line2);
	milieu.add(line3);

	line1.setLayout(layout);
	line2.setLayout(layout);
	line3.setLayout(layout);

	
	//BufferedImage icon = ImageIO.read(getClass().getResourceAsStream("image.png"));
	//ImageIcon icon = new ImageIcon(this.getClass().getResource("image.png"));
	//ImageIcon icon = new ImageIcon("image.png");
	ImageIcon icon = new ImageIcon(this.getClass().getResource(
				"image.png"));
	//JLabel image = new JLabel(new ImageIcon(icon));
	JLabel image = new JLabel(icon);
	haut.add(image);
	
	
	
	
 
	

	log = new JTextArea(17, 27);
	log.setBackground(new Color(102,102,102));
	log.setForeground(Color.white);
	log.setLineWrap(true);
	sbrText = new JScrollPane(log);
	sbrText.setBounds(3, 3, 300, 200);
	//sbrText.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
	

	log.setText("Welcome to LeView! \nCommunication with LeView will appear here...");
	log.setEditable(false);
	//sbrText.add(log);
	haut.add(sbrText);


      JLabel texte = new JLabel("Load a local PDB file:    ");
      line1.add(texte);

	 fichier = new JTextField("",20);
      fichier.setEditable(false);         
      line1.add(fichier);
 
      search = new JButton("  Browse  ");
      line1.add(search);
      search.addActionListener(this);

	JLabel id = new JLabel("Download with PDB id: ");
    line2.add(id);

	iden = new JTextField("",20);      
      line2.add(iden);
 
     

	download = new JButton("Download");
     line2.add(download);
      download.addActionListener(this);
 
 
     open = new JButton(" Run LeView ");
      open.setEnabled(false);
      bas.add(open);
      open.addActionListener(this);
 
     // exit = new JButton("Close");
     // bas.add(exit);
     // exit.addActionListener(this);

	 rb = new JRadioButton("Residue Range (OPTIONAL)", false);
	line3.add(rb);
	rb.addActionListener(this);

	JLabel r1 = new JLabel("Res1:");
      	line3.add(r1);

	res1 = new JTextField("",5);    
	res1.setEditable(false);
      	line3.add(res1);

	JLabel r2 = new JLabel("Res2:");
      	line3.add(r2);

	 res2 = new JTextField("",5);  
	res2.setEditable(false);  
      	line3.add(res2);

	JLabel ch = new JLabel("Chain:");
      	line3.add(ch);

	chain = new JTextField("",5); 
	chain.setEditable(false);   
      	line3.add(chain);

 
	this.add(haut,BorderLayout.NORTH);
	this.add(milieu,BorderLayout.CENTER);
     	this.add(bas,BorderLayout.SOUTH);
 
      this.setVisible(true); 
   }
 
   public void actionPerformed(ActionEvent evenement)
   {
      Object source = evenement.getSource();
 
      if (source == search)
      {
         
         chooser = new JFileChooser();
	chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	FileNameExtensionFilter filter = new FileNameExtensionFilter("PDB files", "pdb", "PDB");
	chooser.setFileFilter(filter);
	String dir=System.getProperty("user.dir");
	File f=new File(dir);
	chooser.setCurrentDirectory(f);
        
         int returnVal = chooser.showOpenDialog(this);
      
         if(returnVal == JFileChooser.APPROVE_OPTION)
         {
            // recuperation du nom du fichier
            fileName = chooser.getSelectedFile().getAbsolutePath();
            // affichage du nom du fichier dans le JTextField
            fichier.setText(fileName);
            // maintenant, on peut cliquer sur le JButton voir
	com=com+"loading "+fileName+"...\n";
	log.setText(com);
            open.setEnabled(true);
         }
 
      }
      else if (source == open)
      {
         
	try {

		
		
		if(!fileName.equals("")){

			if(rb.isSelected()){
				String rr1=res1.getText();
				String rr2=res2.getText();
				String cc=chain.getText();

				if((!rr1.equals(""))&&(!rr2.equals(""))&&(!cc.equals(""))){

				

					try {
					int r1=Integer.parseInt(rr1);
					int r2	=Integer.parseInt(rr2);
					char c=cc.charAt(0);

					if(r1>r2){
						com=com+"ERROR: Problem with residue range. RES1 must be greater than RES2.\n";
						log.setText(com);
						rb.setSelected(false);
						res1.setEditable(false);
						res1.setText(""); 
						res2.setEditable(false); 
						res2.setText("");
						chain.setEditable(false);
						chain.setText("");
					}else{

					File2 f2 = new File2(fileName, r1, r2, c);
					Molecule m = f2.getMolecule();
					
					Traitement2 t = new Traitement2(m);
					LinkedList chains = m.getChains();
					//LinkedList ttt=m.getChains();

					int nb=0;
					for(int i=0;i<chains.size();i++){
						Chain ccc = (Chain) chains.get(i);
						if (ccc.getType().equals("ligand")) {nb++;}
					}

					if(nb==0){
						com=com+"ERROR: Problem with residue range.\n";
						com=com+"RES1="+r1+" RES2="+r2+" ChainID="+c+" does not correspond with a ligand.\n";
						com=com+"Be careful, ChainID is case sensitive!";
						log.setText(com);
					}else{

					
					for (int i = 0; i < chains.size(); i++) {

						Chain ccc = (Chain) chains.get(i);

						//if ((mode == 0) || (mode == 2)) {
							/*if (ccc.getType().equals("ion")) {
								Metal metal = new Metal(ccc, 4.0);
								//new FrameMetal(metal, lm, am, save, colorMode);
								 new MyFrameMetal(metal);
							}*/
						//}
						//if ((mode == 0) || (mode == 1)) {
							if (ccc.getType().equals("ligand")) {
								com=com+"Running LeView  between residues "+r1+" and "+r2+" on chain +"+c+"\n";
								log.setText(com);
								Ligand ligand = new Ligand(ccc, 4.0, 3.3);
								//new Frame(ligand, lm, am, save, hbStyle,displayDist, colorMode);
								 new MyFrame(ligand);
								open.setEnabled(false);
								iden.setText("");
								fichier.setText("");
								rb.setSelected(false);
								res1.setEditable(false);
								res1.setText(""); 
								res2.setEditable(false); 
								res2.setText("");
								chain.setEditable(false);
								chain.setText("");
							}
						//}
					}


					}


					}
					}
					catch (Exception e) {
						com=com+"ERROR: Problem with residue range. Values for RES1 and RES2 must be an integer and a character for Chain.\n";
						log.setText(com);
						rb.setSelected(false);
						res1.setEditable(false);
						res1.setText(""); 
						res2.setEditable(false); 
						res2.setText("");
						chain.setEditable(false);
						chain.setText("");
					
					}
				
				}else{
					com=com+"ERROR: Problem with residue range. You need to give an integer for RES1 and RES2 and a character for Chain.\n";
					log.setText(com);
					rb.setSelected(false);
					res1.setEditable(false);
					res1.setText(""); 
					res2.setEditable(false); 
					res2.setText("");
					chain.setEditable(false);
					chain.setText("");

				}
				
			}
			else{
				com=com+"Running LeView...\n";
				log.setText(com);
				new MainFrame(fileName);
				open.setEnabled(false);
				iden.setText("");
				fichier.setText("");
			}
			
		}
	}
	catch (IOException e) {
 		// Instructions de traitement de l'erreur;
 	}
      }
      else if (source == exit)
      {
         this.dispose();
      }
      else
      {

		if(source==download)
		{
			String ii=iden.getText();
			if(!ii.equals("")){
			try{

			com=com+"Downloading "+ii+" PDB file...\n";
			log.setText(com);

			String tmpo = "http://www.ebi.ac.uk/pdbe-srv/view/files/"+ ii + ".pdb";
			URL monUrl = new URL(tmpo);

			InputStreamReader ins = new InputStreamReader(monUrl.openStream());
			BufferedReader buff = new BufferedReader(ins);
			fileName=ii+".pdb";
			
			String line=buff.readLine();
			if(line.indexOf("HEADER")!=-1){
			PrintWriter out=new PrintWriter(new FileWriter(fileName));
				out.println(line);
				while ((line = buff.readLine()) != null) {
					out.println(line);

				}
			out.close();
			open.setEnabled(true);
			}
			else{
			/* not PDB file ex 4hvv*/
				iden.setText("");
				com=com+"ERROR: "+ii+" is not a correct PDB file\n";
				com=com+"Please check the PDB id on www.ebi.ac.uk/pdbe/\n";
				log.setText(com);
			}

			} catch (IOException ioe) {
			System.out.println("Error --" + ioe.toString());
			iden.setText("");
			com=com+"ERROR: PDB id +"+ ii +"is not found\n";
			com=com+"Please check the PDB id on www.ebi.ac.uk/pdbe/\n";
			log.setText(com);
			}
			}

		}
		else{

			if(source==rb){
			if(rb.isSelected()){
				res1.setEditable(true); 
				res2.setEditable(true); 
				chain.setEditable(true); 
			}else{
				res1.setEditable(false);
				res1.setText(""); 
				res2.setEditable(false); 
				res2.setText("");
				chain.setEditable(false);
				chain.setText("");
			}
			}
			else{

			}

		}
      }
   }
 
  
 
   /*
    * point d'entre du programme
    */
   public static void main(String[] arguments)
   {
      Chooser exemple = new Chooser();
   }

 
   
}
