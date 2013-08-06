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
 * Upper symmetrix dense matrix. It has the same storage layout as the
 * {@link no.uib.cipr.matrix.DenseMatrix DenseMatrix}, but only refers to
 * elements above or on the main diagonal. The remaining elements are never
 * accessed nor changed, and is known only by symmetry.
 */
public class UpperSymmDenseMatrix extends AbstractSymmDenseMatrix {

    /**
     * Constructor for UpperSymmDenseMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public UpperSymmDenseMatrix(int n) {
        super(n, UpLo.Upper);
    }

    /**
     * Constructor for UpperSymmDenseMatrix
     * 
     * @param A
     *            Matrix to copy. It must be a square matrix, and only the upper
     *            triangular part is copied
     */
    public UpperSymmDenseMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for UpperSymmDenseMatrix
     * 
     * @param A
     *            Matrix to copy. It must be a square matrix, and only the upper
     *            triangular part is copied
     * @param deep
     *            If false, a shallow copy is made. In that case, <code>A</code>
     *            must be a dense matrix
     */
    public UpperSymmDenseMatrix(Matrix A, boolean deep) {
        super(A, deep, UpLo.Upper);
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
    public UpperSymmDenseMatrix copy() {
        return new UpperSymmDenseMatrix(this);
    }

    @Override
    void copy(Matrix A) {
        for (MatrixEntry e : A)
            if (e.row() <= e.column())
                set(e.row(), e.column(), e.get());
    }

}
