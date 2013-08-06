/*
 * Copyright (C) 2003-2006 Bjørn-Ove Heimsund
 * 
 * This file is part of MTJ.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
package no.uib.cipr.matrix;

/**
 * The job the singular value solvers are to do. This only limits which singular vectors
 * are computed, all the singular values are always computed
 */
enum JobSVD {
	/** Compute all of the singular vectors */
	All,

	/** Do not compute any singular vectors */
	None,

	/**
	 * Overwrite passed data. For an <code>M*N</code> matrix, this either overwrites the
	 * passed matrix with as many singular vectors as there is room for. Details depend on
	 * the actual algorithm
	 */
	Overwrite,

	/**
	 * Compute parts of the singular vectors. For an <code>M*N</code> matrix, this
	 * computes <code>min(M,N)</code> singular vectors
	 */
	Part;

	/**
	 * @return the netlib character version of this designation, for use with F2J.
	 */
	public String netlib() {
		switch (this) {
		case All:
			return "A";
		case Part:
			return "S";
		case Overwrite:
			return "O";
		default:
			return "N";
		}
	}

}