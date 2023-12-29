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
 This class defines node for the Graph class.
 */

import java.util.*;

public class Node {

	public int n1;
	public int n2;

	public Node() {
	}

	public Node(int u, int v) {
		n1 = u;
		n2 = v;
	}

	public int getOne() {
		return n1;
	}

	public int getTwo() {
		return n2;
	}

	public void print() {
		System.out.println("le premier " + getOne() + "le second " + getTwo());
	}

}
