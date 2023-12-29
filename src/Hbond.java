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
 This class defines an hydrogen bonds with a donnor and an acceptor.
 */


import java.io.*;
import java.util.*;

public class Hbond {

	Atom don;
	Atom acc;
	Residue res;
	int type;
	boolean flag;
	LinkedList linked;
	double distance;
	boolean selected;

	public Hbond() {
	}

	public Hbond(Atom a1, Atom a2, Residue r) {
		don = a1;
		acc = a2;
		res = r;
		flag = false;
		linked = new LinkedList();
		distance = 0;
		selected = false;
	}

	public Hbond(Atom a1, Atom a2, Residue r, double distance) {
		don = a1;
		acc = a2;
		res = r;
		flag = false;
		linked = new LinkedList();
		this.distance = distance;
		selected = false;
	}

	public void addLinked(Hbond h) {
		linked.add(h);
	}

	public void setLinked(LinkedList l) {
		linked = l;
	}

	public boolean getSelected() {
		return selected;
	}

	public void setSelected(boolean b) {
		selected = b;
	}

	public double getDistance() {
		return distance;
	}

	public LinkedList getLinked() {
		return linked;
	}

	public Atom getDon() {
		return don;
	}

	public boolean getFlag() {
		return flag;
	}

	public void setFlag(boolean b) {
		flag = b;
	}

	public Atom getAcc() {
		return acc;
	}

	public Residue getParent() {
		return res;
	}

	public void setOrder(Residue r) {
		res = r;
	}

	public Atom getAtomRes() {
		if (don.getParent() == res)
			return don;
		else
			return acc;
	}

	public Atom getAtomLig() {
		if (don.getParent() == res)
			return acc;
		else
			return don;
	}

	public void rotation(double angle) {
		
		/* angle in degrees */

		angle = Math.toRadians(angle);

		Atom at = getAtomRes();
		Atom centre = getAtomLig();

		double OAx = at.getFx() - centre.getFx();
		double OAy = at.getFy() - centre.getFy();

		double resx = centre.getFx() + Math.cos(angle) * OAx - Math.sin(angle)
				* OAy;
		double resy = centre.getFy() + Math.sin(angle) * OAx + Math.cos(angle)
				* OAy;

		at.setFx(resx);
		at.setFy(resy);

	}

	public boolean equals(Hbond b) {
		if ((b.getDon().getId() == don.getId())
				&& (b.getAcc().getId() == acc.getId()))
			return true;
		if ((b.getDon().getId() == acc.getId())
				&& (b.getAcc().getId() == don.getId()))
			return true;
		return false;
	}

	public void print() {
		System.out.println("HBOND " + don.getId() + " " + acc.getId() + " "
				+ res.getName());
	}

}
