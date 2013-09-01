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

import com.github.fommil.netlib.BLAS;
import com.github.fommil.netlib.LAPACK;
import org.netlib.util.intW;

import java.util.Iterator;


/**
 * Partial implementation of a triangular, packed matrix
 */
abstract class AbstractTriangPackMatrix extends AbstractPackMatrix {

    /**
     * Upper or lower triangular
     */
    UpLo uplo;

    /**
     * If the matrix is unit diagonal or not unit
     */
    Diag diag;

    /**
     * Constructor for AbstractTriangPackMatrix
     */
    AbstractTriangPackMatrix(int n, UpLo uplo, Diag diag) {
        super(n);
        this.uplo = uplo;
        this.diag = diag;
    }

    /**
     * Constructor for AbstractTriangPackMatrix
     */
    AbstractTriangPackMatrix(Matrix A, UpLo uplo, Diag diag) {
        this(A, false, uplo, diag);
    }

    /**
     * Constructor for AbstractTriangPackMatrix
     */
    AbstractTriangPackMatrix(Matrix A, boolean deep, UpLo uplo, Diag diag) {
        super(A, deep);
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
        BLAS.getInstance().dtpmv(uplo.netlib(), Transpose.NoTranspose.netlib(), diag.netlib(),
        	numRows, data, yd, 1);

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

        // y = A*z
        BLAS.getInstance().dtpmv(uplo.netlib(), Transpose.Transpose.netlib(), diag.netlib(),
        	numRows, data, yd, 1);

        return y;
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

        checkSolve(B, X);

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        intW info = new intW(0);
        LAPACK.getInstance().dtptrs(uplo.netlib(), trans.netlib(), diag.netlib(), numRows,
                X.numColumns(), data, Xd, Matrices.ld(numRows), info);

        if (info.val > 0)
            throw new MatrixSingularException();
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return X;
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new TriangPackMatrixIterator();
    }

    private class TriangPackMatrixIterator extends RefMatrixIterator {

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
