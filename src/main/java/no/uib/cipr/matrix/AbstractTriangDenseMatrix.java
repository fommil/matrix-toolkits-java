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

import java.util.Iterator;

import com.github.fommil.netlib.BLAS;
import com.github.fommil.netlib.LAPACK;
import org.netlib.util.intW;


/**
 * Partial implementation of a triangular, dense matrix
 */
abstract class AbstractTriangDenseMatrix extends AbstractDenseMatrix {

    /**
     * If the matrix is upper triangular
     */
    UpLo uplo;

    /**
     * If the matrix is unit diagonal or not unit
     */
    Diag diag;

    /**
     * Leading dimension of the matrix
     */
    int ld;

    /**
     * Constructor for AbstractTriangDenseMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    AbstractTriangDenseMatrix(int n, UpLo uplo, Diag diag) {
        super(n, n);
        ld = n;
        this.uplo = uplo;
        this.diag = diag;
    }

    /**
     * Constructor for AbstractTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from
     */
    AbstractTriangDenseMatrix(Matrix A, UpLo uplo, Diag diag) {
        this(A, Math.min(A.numRows(), A.numColumns()), uplo, diag);
    }

    /**
     * Constructor for AbstractTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from
     * @param deep
     *            If true, <code>A</code> is copied, else a shallow copy is
     *            made and the matrices share underlying storage. For this,
     *            <code>A</code> must be a dense matrix
     */
    AbstractTriangDenseMatrix(Matrix A, boolean deep, UpLo uplo, Diag diag) {
        this(A, Math.min(A.numRows(), A.numColumns()), deep, uplo, diag);
    }

    /**
     * Constructor for AbstractTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from
     * @param k
     *            Size of matrix to refer.
     *            <code>k&lt;min(numRows,numColumns)</code>
     */
    AbstractTriangDenseMatrix(Matrix A, int k, UpLo uplo, Diag diag) {
        this(A, k, true, uplo, diag);
    }

    /**
     * Constructor for AbstractTriangDenseMatrix
     * 
     * @param A
     *            Matrix to copy from
     * @param k
     *            Size of matrix to refer.
     *            <code>k&lt;min(numRows,numColumns)</code>
     * @param deep
     *            If true, <code>A</code> is copied, else a shallow copy is
     *            made and the matrices share underlying storage. For this,
     *            <code>A</code> must be a dense matrix
     */
    AbstractTriangDenseMatrix(Matrix A, int k, boolean deep, UpLo uplo,
            Diag diag) {
        super(A, deep);

        if (deep && !A.isSquare())
            throw new IllegalArgumentException("deep && !A.isSquare()");

        ld = A.numRows();
        numRows = numColumns = k;
        this.uplo = uplo;
        this.diag = diag;
    }

    @Override
    public Vector mult(double alpha, Vector x, Vector y) {
        if (!(y instanceof DenseVector))
            return super.mult(alpha, x, y);

        checkMultAdd(x, y);

        double[] yd = ((DenseVector) y).getData();

        // y = alpha*x
        y.set(alpha, x);

        // y = A*z
        BLAS.getInstance().dtrmv(uplo.netlib(), Transpose.NoTranspose.netlib(), diag.netlib(),
        	numRows, data, Math.max(1, ld), yd, 1);

        return y;
    }

    @Override
    public Vector transMult(double alpha, Vector x, Vector y) {
        if (!(y instanceof DenseVector))
            return super.transMult(alpha, x, y);

        checkTransMultAdd(x, y);

        double[] yd = ((DenseVector) y).getData();

        // y = alpha*x
        y.set(alpha, x);

        // y = A'*y
        BLAS.getInstance().dtrmv(uplo.netlib(), Transpose.Transpose.netlib(), diag.netlib(),
        	numRows, data, Math.max(1, ld), yd, 1);

        return y;
    }

    @Override
    public Matrix mult(double alpha, Matrix B, Matrix C) {
        if (!(C instanceof DenseMatrix))
            return super.mult(alpha, B, C);

        checkMultAdd(B, C);

        double[] Cd = ((DenseMatrix) C).getData();

        C.set(B);

        // C = alpha*A*C
        BLAS.getInstance().dtrmm(Side.Left.netlib(), uplo.netlib(), Transpose.NoTranspose.netlib(),
        	diag.netlib(), C.numRows(), C.numColumns(), alpha, data, Math.max(1, ld), Cd,
        	Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public Matrix transAmult(double alpha, Matrix B, Matrix C) {
        if (!(C instanceof DenseMatrix))
            return super.transAmult(alpha, B, C);

        checkTransAmultAdd(B, C);

        double[] Cd = ((DenseMatrix) C).getData();

        C.set(B);

        // C = alpha*A'*C
        BLAS.getInstance().dtrmm(Side.Left.netlib(), uplo.netlib(), Transpose.Transpose.netlib(),
        	diag.netlib(), C.numRows(), C.numColumns(), alpha, data, Math.max(1, ld), Cd,
        	Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        return solve(B, X, Transpose.NoTranspose);
    }

    @Override
    public Vector solve(Vector b, Vector x) {
        DenseMatrix B = new DenseMatrix(b, false), X = new DenseMatrix(x, false);
        solve(B, X);
        return x;
    }

    @Override
    public Matrix transSolve(Matrix B, Matrix X) {
        return solve(B, X, Transpose.Transpose);
    }

    @Override
    public Vector transSolve(Vector b, Vector x) {
        DenseMatrix B = new DenseMatrix(b, false), X = new DenseMatrix(x, false);
        transSolve(B, X);
        return x;
    }

    Matrix solve(Matrix B, Matrix X, Transpose trans) {
        if (!(X instanceof DenseMatrix))
            throw new UnsupportedOperationException("X must be a DenseMatrix");

        // Different argument checking to support Hessenberg type matrices for
        // solvers such as GMRES
        if (B.numRows() < numRows)
            throw new IllegalArgumentException("B.numRows() < A.numRows()");
        if (B.numColumns() != X.numColumns())
            throw new IllegalArgumentException(
                    "B.numColumns() != X.numColumns()");
        if (X.numRows() < numRows)
            throw new IllegalArgumentException("X.numRows() < A.numRows()");

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        intW info = new intW(0);
        LAPACK.getInstance().dtrtrs(uplo.netlib(), trans.netlib(), diag.netlib(), numRows,
                X.numColumns(), data, Math.max(1, ld), Xd, Matrices.ld(numRows), info);

        if (info.val > 0)
            throw new MatrixSingularException();
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return X;
    }

    @Override
    int getIndex(int row, int column) {
        check(row, column);
        return row + column * Math.max(ld, numRows);
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new TriangDenseMatrixIterator();
    }

    private class TriangDenseMatrixIterator extends RefMatrixIterator {

        @Override
        public MatrixEntry next() {
            entry.update(row, column);

            if (uplo == UpLo.Lower)
                if (row < numRows - 1)
                    row++;
                else {
                    column++;
                    row = column;
                }
            else { // uplo == UpLo.Upper
                if (row < column)
                    row++;
                else {
                    column++;
                    row = 0;
                }
            }

            return entry;
        }

    }

}
