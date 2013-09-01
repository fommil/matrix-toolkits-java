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

/**
 * Partial implementation of a triangular, banded matrix
 */
abstract class AbstractTriangBandMatrix extends AbstractBandMatrix {

    /**
     * Upper or lower part stored?
     */
    private UpLo uplo;

    /**
     * Diagonal stored or not?
     */
    private Diag diag;

    /**
     * Diagonals in relevant band
     */
    int kd;

    /**
     * Constructor for AbstractTriangBandMatrix
     */
    AbstractTriangBandMatrix(int n, int kl, int ku, UpLo uplo, Diag diag) {
        super(n, kl, ku);
        kd = Math.max(kl, ku);
        this.uplo = uplo;
        this.diag = diag;
    }

    /**
     * Constructor for AbstractTriangBandMatrix
     */
    AbstractTriangBandMatrix(Matrix A, int kl, int ku, UpLo uplo, Diag diag) {
        this(A, kl, ku, true, uplo, diag);
    }

    /**
     * Constructor for AbstractTriangBandMatrix
     */
    AbstractTriangBandMatrix(Matrix A, int kl, int ku, boolean deep, UpLo uplo,
            Diag diag) {
        super(A, kl, ku, deep);
        kd = Math.max(kl, ku);
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

        // y = A*y
        BLAS.getInstance().dtbmv(uplo.netlib(), Transpose.NoTranspose.netlib(), diag.netlib(),
        	numRows, kd, data, kd + 1, yd, 1);

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

        // y = A*y
        BLAS.getInstance().dtbmv(uplo.netlib(), Transpose.Transpose.netlib(), diag.netlib(),
        	numRows, kd, data, kd + 1, yd, 1);

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
        LAPACK.getInstance().dtbtrs(uplo.netlib(), trans.netlib(), diag.netlib(),
        	numRows, kd, X.numColumns(), data, Matrices.ld(kd + 1), Xd, Matrices.ld(n),
        	info);

        if (info.val > 0)
            throw new MatrixSingularException();
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return X;
    }

}
