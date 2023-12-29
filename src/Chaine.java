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
 This class defines a ligand chain. A chain is composed at least of 3 non-terminal heavy atoms.
 */

import java.io.*;
import java.util.*;

public class Chaine {

	LinkedList atomId;

	public Chaine() {
	}

	public Chaine(LinkedList l) {
		atomId = l;
	}

	public LinkedList getAtoms() {
		return atomId;
	}

	public int getSize() {
		return atomId.size();
	}

	public void print() {

		System.out.println("CHAINE " + atomId);
	}

}
