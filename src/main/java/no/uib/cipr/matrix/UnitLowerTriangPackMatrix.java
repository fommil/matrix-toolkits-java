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
 * Unit lower triangular packed matrix. Same storage as
 * {@link no.uib.cipr.matrix.LowerTriangPackMatrix LowerTriangPackMatrix}, but
 * the main diagonal is assumed to be all ones.
 */
public class UnitLowerTriangPackMatrix extends LowerTriangPackMatrix {

    /**
     * Constructor for UnitLowerTriangPackMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public UnitLowerTriangPackMatrix(int n) {
        super(n, Diag.Unit);
    }

    /**
     * Constructor for UnitLowerTriangPackMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the entries of the relevant
     *            part are copied
     */
    public UnitLowerTriangPackMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for UnitLowerTriangPackMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the entries of the relevant
     *            part are copied
     * @param deep
     *            True if the copy is deep, else false (giving a shallow copy).
     *            For shallow copies, <code>A</code> must be a packed matrix
     */
    public UnitLowerTriangPackMatrix(Matrix A, boolean deep) {
        super(A, deep, Diag.Unit);
    }

    @Override
    public void add(int row, int column, double value) {
        if (column == row)
            throw new IllegalArgumentException("column == row");
        super.add(row, column, value);
    }

    @Override
    public double get(int row, int column) {
        if (column == row)
            return 1;
        return super.get(row, column);
    }

    @Override
    public void set(int row, int column, double value) {
        if (column == row)
            throw new IllegalArgumentException("column == row");
        super.set(row, column, value);
    }

    @Override
    void copy(Matrix A) {
        for (MatrixEntry e : A)
            if (e.row() > e.column())
                set(e.row(), e.column(), e.get());
    }

    @Override
    public UnitLowerTriangPackMatrix copy() {
        return new UnitLowerTriangPackMatrix(this);
    }

}
