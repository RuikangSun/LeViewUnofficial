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
 This class defines the applet.
 */


import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JApplet;
import javax.swing.JButton;


import java.io.*;
import java.util.*;


public class MyApplet extends JApplet {
	
	private static final long serialVersionUID = 1L;
	private JButton boutonDemarrerArreter;

	
	public MyApplet() {
		super();
	}

	
	public void init()   {
		super.init();
		
		setBackground(Color.white);
		add (boutonDemarrerArreter = new JButton ("Click"));
		boutonDemarrerArreter.setName("Demarrer");

		boutonDemarrerArreter.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				if (((JButton) e.getSource()).getName() == "Demarrer")
				{

				try{
				MainFrame frame=new MainFrame(getParameter("file"));
				}catch(java.io.IOException exp){ exp.printStackTrace();}
					
				}
			}
		});
		
	}
}

