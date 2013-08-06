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
 * Lower triangular packed matrix. In contrast with
 * {@link no.uib.cipr.matrix.LowerTriangDenseMatrix LowerTriangDenseMatrix},
 * this matrix exploits the sparsity by only storing about half the matrix. As
 * such, the triangular matrix
 * <p>
 * <table border="1">
 * <tr>
 * <td>a<sub>11</sub></td>
 * <td>&nbsp;</td>
 * <td>&nbsp;</td>
 * <td>&nbsp;</td>
 * </tr>
 * <tr>
 * <td>a<sub>21</sub></td>
 * <td>a<sub>22</sub></td>
 * <td>&nbsp;</td>
 * <td>&nbsp;</td>
 * </tr>
 * <tr>
 * <td>a<sub>31</sub></td>
 * <td>a<sub>32</sub></td>
 * <td>a<sub>33</sub></td>
 * <td>&nbsp;</td>
 * </tr>
 * <tr>
 * <td>a<sub>41</sub></td>
 * <td>a<sub>42</sub></td>
 * <td>a<sub>43</sub></td>
 * <td>a<sub>44</sub></td>
 * </tr>
 * </table>
 * </p>
 * <p>
 * is packed as follows:
 * </p>
 * <p>
 * <table border="1">
 * <tr>
 * <td>a<sub>11</sub></td>
 * <td>a<sub>21</sub></td>
 * <td>a<sub>31</sub></td>
 * <td>a<sub>41</sub></td>
 * <td>a<sub>22</sub></td>
 * <td>a<sub>32</sub></td>
 * <td>a<sub>42</sub></td>
 * <td>a<sub>33</sub></td>
 * <td>a<sub>43</sub></td>
 * <td>a<sub>44</sub></td>
 * </tr>
 * </table>
 * </p>
 */
public class LowerTriangPackMatrix extends AbstractTriangPackMatrix {

    /**
     * Constructor for LowerTriangPackMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public LowerTriangPackMatrix(int n) {
        super(n, UpLo.Lower, Diag.NonUnit);
    }

    /**
     * Constructor for LowerTriangPackMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    LowerTriangPackMatrix(int n, Diag diag) {
        super(n, UpLo.Lower, diag);
    }

    /**
     * Constructor for LowerTriangPackMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the entries of the relevant
     *            part are copied
     */
    public LowerTriangPackMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for LowerTriangPackMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the entries of the relevant
     *            part are copied
     * @param deep
     *            True if the copy is deep, else false (giving a shallow copy).
     *            For shallow copies, <code>A</code> must be a packed matrix
     */
    public LowerTriangPackMatrix(Matrix A, boolean deep) {
        super(A, deep, UpLo.Lower, Diag.NonUnit);
    }

    /**
     * Constructor for LowerTriangPackMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the entries of the relevant
     *            part are copied
     * @param deep
     *            True if the copy is deep, else false (giving a shallow copy).
     *            For shallow copies, <code>A</code> must be a packed matrix
     */
    LowerTriangPackMatrix(Matrix A, boolean deep, Diag diag) {
        super(A, deep, UpLo.Lower, diag);
    }

    @Override
    public void add(int row, int column, double value) {
        if (column > row)
            throw new IllegalArgumentException("column > row");
        data[getIndex(row, column)] += value;
    }

    @Override
    public void set(int row, int column, double value) {
        if (column > row)
            throw new IllegalArgumentException("column > row");
        data[getIndex(row, column)] = value;
    }

    @Override
    public double get(int row, int column) {
        if (column > row)
            return 0;
        return data[getIndex(row, column)];
    }

    /**
     * Checks the row and column indices, and returns the linear data index
     */
    int getIndex(int row, int column) {
        check(row, column);
        return row + (2 * n - (column + 1)) * column / 2;
    }

    @Override
    void copy(Matrix A) {
        for (MatrixEntry e : A)
            if (e.row() >= e.column())
                set(e.row(), e.column(), e.get());
    }

    @Override
    public LowerTriangPackMatrix copy() {
        return new LowerTriangPackMatrix(this);
    }

}
