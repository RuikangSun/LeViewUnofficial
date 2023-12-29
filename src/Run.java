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
 This class is the main class.
 */

import java.io.*;
import java.util.*;

public class Run {

	public static void main(String[] args) throws IOException

	{
		
		if(args.length!=0){
			String fic=args[0];
		
			if (fic.endsWith(".pdb")) {
				MainFrame mf=new MainFrame(fic);

			}
			else {
			System.out.println("no pdb file found");
			System.exit(0);
		}
		}
		else {
			System.out.println("no pdb file found");
			System.exit(0);
		}

	}
}
