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
 * Upper symmetrical banded matrix. The same storage as
 * {@link no.uib.cipr.matrix.BandMatrix BandMatrix}, but without subdiagonals.
 * Lower part of the matrix is implictly known by symmetry
 */
public class UpperSymmBandMatrix extends AbstractSymmBandMatrix {

    /**
     * Constructor for UpperSymmBandMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     * @param kd
     *            Number of bands off the main diagonal (off diagonals)
     */
    public UpperSymmBandMatrix(int n, int kd) {
        super(n, 0, kd, UpLo.Upper);
    }

    /**
     * Constructor for UpperSymmBandMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the parts of <code>A</code>
     *            that lie within the allocated band are copied over, the rest
     *            is ignored
     * @param kd
     *            Number of bands off the main diagonal (off diagonals)
     */
    public UpperSymmBandMatrix(Matrix A, int kd) {
        this(A, kd, true);
    }

    /**
     * Constructor for UpperSymmBandMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the parts of <code>A</code>
     *            that lie within the allocated band are copied over, the rest
     *            is ignored
     * @param kd
     *            Number of bands off the main diagonal (off diagonals)
     * @param deep
     *            True for a deep copy. For shallow copies, <code>A</code>
     *            must be a banded matrix
     */
    public UpperSymmBandMatrix(Matrix A, int kd, boolean deep) {
        super(A, 0, kd, deep, UpLo.Upper);
    }

    @Override
    public void add(int row, int column, double value) {
        if (row <= column)
            super.add(row, column, value);
    }

    @Override
    public double get(int row, int column) {
        if (row > column)
            return super.get(column, row);
        return super.get(row, column);
    }

    @Override
    public void set(int row, int column, double value) {
        if (row <= column)
            super.set(row, column, value);
    }

    @Override
    public UpperSymmBandMatrix copy() {
        return new UpperSymmBandMatrix(this, kd);
    }

}
