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
 This class defines a chemical bond by 2 atoms and order.
 */

import java.io.*;
import java.util.*;

public class Bond {

	Atom a1;
	Atom a2;
	int order;

	public Bond() {
	}

	public Bond(Atom a1, Atom a2, int order) {
		this.a1 = a1;
		this.a2 = a2;
		this.order = order;
	}

	public Atom getFirst() {
		return a1;
	}

	public Atom getSecond() {
		return a2;
	}

	public int getOrder() {
		return order;
	}

	public void setOrder(int order) {
		this.order = order;
	}

	public boolean equals(Bond b) {
		if ((b.getFirst().getId() == a1.getId())
				&& (b.getSecond().getId() == a2.getId()))
			return true;
		if ((b.getFirst().getId() == a2.getId())
				&& (b.getSecond().getId() == a1.getId()))
			return true;
		return false;
	}

	public void print() {
		System.out.println("BOND " + a1.getId() + " " + a2.getId() + " "
				+ order);
	}

}
