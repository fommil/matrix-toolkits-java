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
 * Lower triangular dense matrix. It has the same storage layout as the
 * {@link no.uib.cipr.matrix.DenseMatrix DenseMatrix}, but only refers to
 * elements below or on the main diagonal. The remaining elements are assumed to
 * be zero, but since they are never accessed, they need not be.
 */
public class LowerTriangDenseMatrix extends AbstractTriangDenseMatrix {

    /**
     * Constructor for LowerTriangDenseMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public LowerTriangDenseMatrix(int n) {
        super(n, UpLo.Lower, Diag.NonUnit);
    }

    /**
     * Constructor for LowerTriangDenseMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    LowerTriangDenseMatrix(int n, Diag diag) {
        super(n, UpLo.Lower, diag);
    }

    /**
     * Constructor for LowerTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from. Only the lower triangular part is copied
     */
    public LowerTriangDenseMatrix(Matrix A) {
        this(A, Math.min(A.numRows(), A.numColumns()));
    }

    /**
     * Constructor for LowerTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from. Only the lower triangular part is copied
     * @param deep
     *            If true, <code>A</code> is copied, else a shallow copy is
     *            made and the matrices share underlying storage. For this,
     *            <code>A</code> must be a dense matrix
     */
    public LowerTriangDenseMatrix(Matrix A, boolean deep) {
        this(A, Math.min(A.numRows(), A.numColumns()), deep);
    }

    /**
     * Constructor for LowerTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from. Only the lower triangular part is copied
     * @param deep
     *            If true, <code>A</code> is copied, else a shallow copy is
     *            made and the matrices share underlying storage. For this,
     *            <code>A</code> must be a dense matrix
     */
    LowerTriangDenseMatrix(Matrix A, boolean deep, Diag diag) {
        this(A, Math.min(A.numRows(), A.numColumns()), deep, diag);
    }

    /**
     * Constructor for LowerTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from. Only the lower triangular part is copied
     * @param k
     *            Size of matrix to refer.
     *            <code>k&lt;min(numRows,numColumns)</code>
     */
    public LowerTriangDenseMatrix(Matrix A, int k) {
        this(A, k, true);
    }

    /**
     * Constructor for LowerTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from. Only the lower triangular part is copied
     * @param k
     *            Size of matrix to refer.
     *            <code>k&lt;min(numRows,numColumns)</code>
     * @param deep
     *            If true, <code>A</code> is copied, else a shallow copy is
     *            made and the matrices share underlying storage. For this,
     *            <code>A</code> must be a dense matrix
     */
    public LowerTriangDenseMatrix(Matrix A, int k, boolean deep) {
        super(A, k, deep, UpLo.Lower, Diag.NonUnit);
    }

    /**
     * Constructor for LowerTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from. Only the lower triangular part is copied
     * @param k
     *            Size of matrix to refer.
     *            <code>k&lt;min(numRows,numColumns)</code>
     * @param deep
     *            If true, <code>A</code> is copied, else a shallow copy is
     *            made and the matrices share underlying storage. For this,
     *            <code>A</code> must be a dense matrix
     */
    LowerTriangDenseMatrix(Matrix A, int k, boolean deep, Diag diag) {
        super(A, k, deep, UpLo.Lower, diag);
    }

    @Override
    public void add(int row, int column, double value) {
        if (column > row)
            throw new IllegalArgumentException("column > row");
        super.add(row, column, value);
    }

    @Override
    public double get(int row, int column) {
        if (column > row)
            return 0;
        return super.get(row, column);
    }

    @Override
    public void set(int row, int column, double value) {
        if (column > row)
            throw new IllegalArgumentException("column > row");
        super.set(row, column, value);
    }

    @Override
    public LowerTriangDenseMatrix copy() {
        return new LowerTriangDenseMatrix(this);
    }

    @Override
    void copy(Matrix A) {
        for (MatrixEntry e : A)
            if (e.row() >= e.column())
                set(e.row(), e.column(), e.get());
    }

  @Override
  public Matrix set(Matrix A) {
    zero();
    copy(A);
    return this;
  }

}
