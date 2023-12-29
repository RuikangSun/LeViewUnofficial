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
 This class defines a secondary structure.
 */

import java.io.*;
import java.util.*;

public class Structure {

	char type;
	char cbeg;
	char cend;
	int beg;
	int end;
	int id;

	public Structure() {

	}

	public int getBeg() {
		return beg;
	}

	public int getEnd() {
		return end;
	}

	public char getCbeg() {
		return cbeg;
	}

	public char getCend() {
		return cend;
	}

	public char getType() {
		return type;
	}

	public void setType(char c) {
		type = c;
	}

	public void setId(int c) {
		id = c;
	}

	public void setBeg(int c) {
		beg = c;
	}

	public void setEnd(int c) {
		end = c;
	}

	public void setCbeg(char c) {
		cbeg = c;
	}

	public void setCend(char c) {
		cend = c;
	}

	public boolean containRes(Residue res) {
		Chain c = res.getParent();
		if (c.getName() != cbeg)
			return false;

		if ((beg <= res.getIdInt()) && res.getIdInt() <= end)
			return true;
		return false;
	}

	public char containRes2(Residue res) {
		Chain c = res.getParent();
		if (c.getName() != cbeg)
			return 'K';

		if ((beg <= res.getIdInt()) && res.getIdInt() <= end)
			return type;
		return 'K';
	}

	public void print() {
		System.out.println("STRUCTURE " + type + " " + cbeg + " " + beg + " "
				+ cend + " " + end);
	}

}
