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
 * Upper symmetric packed matrix. Same storage as
 * {@link no.uib.cipr.matrix.UpperTriangPackMatrix UpperTriangPackMatrix}, but
 * the lower triangular part is known by symmetry.
 */
public class UpperSymmPackMatrix extends AbstractSymmPackMatrix {

    /**
     * Constructor for UpperSymmPackMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public UpperSymmPackMatrix(int n) {
        super(n, UpLo.Upper);
    }

    /**
     * Constructor for UpperSymmPackMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the entries of the relevant
     *            part are copied
     */
    public UpperSymmPackMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for UpperSymmPackMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the entries of the relevant
     *            part are copied
     * @param deep
     *            True if the copy is deep, else false (giving a shallow copy).
     *            For shallow copies, <code>A</code> must be a packed matrix
     */
    public UpperSymmPackMatrix(Matrix A, boolean deep) {
        super(A, deep, UpLo.Upper);
    }

    @Override
    public void add(int row, int column, double value) {
        if (row <= column)
            data[getIndex(row, column)] += value;
    }

    @Override
    public void set(int row, int column, double value) {
        if (row <= column)
            data[getIndex(row, column)] = value;
    }

    @Override
    public double get(int row, int column) {
        if (row <= column)
            return data[getIndex(row, column)];
        return data[getIndex(column, row)];
    }

    /**
     * Checks the row and column indices, and returns the linear data index
     */
    int getIndex(int row, int column) {
        check(row, column);
        return row + (column + 1) * column / 2;
    }

    @Override
    void copy(Matrix A) {
        for (MatrixEntry e : A)
            if (e.row() <= e.column())
                set(e.row(), e.column(), e.get());
    }

    @Override
    public UpperSymmPackMatrix copy() {
        return new UpperSymmPackMatrix(this);
    }

}
