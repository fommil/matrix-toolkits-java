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
 * Unit lower triangular banded matrix. The same storage as
 * {@link no.uib.cipr.matrix.LowerTriangBandMatrix LowerTriangBandMatrix}, but
 * the main diagonal is assumed to be all ones. It is still allocated, however.
 */
public class UnitLowerTriangBandMatrix extends LowerTriangBandMatrix {

    /**
     * Constructor for UnitLowerTriangBandMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     * @param kd
     *            Number of bands below the main diagonal (subdiagonals)
     */
    public UnitLowerTriangBandMatrix(int n, int kd) {
        super(n, kd, Diag.Unit);
    }

    /**
     * Constructor for UnitLowerTriangBandMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the parts of <code>A</code>
     *            that lie within the allocated band are copied over, the rest
     *            is ignored (including main diagonal entries)
     * @param kd
     *            Number of bands below the main diagonal (subdiagonals)
     */
    public UnitLowerTriangBandMatrix(Matrix A, int kd) {
        this(A, kd, true);
    }

    /**
     * Constructor for UnitLowerTriangBandMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the parts of <code>A</code>
     *            that lie within the allocated band are copied over, the rest
     *            is ignored (including main diagonal entries)
     * @param kd
     *            Number of bands below the main diagonal (subdiagonals)
     * @param deep
     *            True for a deep copy. For shallow copies, <code>A</code>
     *            must be a banded matrix
     */
    public UnitLowerTriangBandMatrix(Matrix A, int kd, boolean deep) {
        super(A, kd, deep, Diag.Unit);
    }

    @Override
    public void add(int row, int column, double value) {
        if (row == column)
            throw new IllegalArgumentException("row == column");
        super.add(row, column, value);
    }

    @Override
    public double get(int row, int column) {
        if (row == column)
            return 1;
        return super.get(row, column);
    }

    @Override
    public void set(int row, int column, double value) {
        if (row == column)
            throw new IllegalArgumentException("row == column");
        super.set(row, column, value);
    }

    @Override
    public UnitLowerTriangBandMatrix copy() {
        return new UnitLowerTriangBandMatrix(this, kl);
    }

    @Override
    public Matrix zero() {
        throw new UnsupportedOperationException();
    }

    @Override
    void copy(Matrix A) {
        for (MatrixEntry e : A)
            if (inBand(e.row(), e.column()) && e.row() != e.column())
                set(e.row(), e.column(), e.get());
    }

}
