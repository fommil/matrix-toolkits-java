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
 * Partial implementation of a symmetrical, banded matrix
 */
abstract class AbstractSymmBandMatrix extends AbstractBandMatrix {

    /**
     * Upper or lower part stored?
     */
    private UpLo uplo;

    /**
     * Diagonals in relevant band
     */
    int kd;

    /**
     * Constructor for AbstractSymmBandMatrix
     */
    AbstractSymmBandMatrix(int n, int kl, int ku, UpLo uplo) {
        super(n, kl, ku);
        kd = Math.max(kl, ku);
        this.uplo = uplo;
    }

    /**
     * Constructor for AbstractSymmBandMatrix
     */
    AbstractSymmBandMatrix(Matrix A, int kl, int ku, UpLo uplo) {
        this(A, kl, ku, true, uplo);
    }

    /**
     * Constructor for AbstractSymmBandMatrix
     */
    AbstractSymmBandMatrix(Matrix A, int kl, int ku, boolean deep, UpLo uplo) {
        super(A, kl, ku, deep);
        kd = Math.max(kl, ku);
        this.uplo = uplo;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.multAdd(alpha, x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dsbmv(uplo.netlib(), numRows, kd, alpha, data, kd + 1, xd, 1, 1, yd, 1);

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        return multAdd(alpha, x, y);
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new BandMatrixIterator(kd, kd);
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        if (!(X instanceof DenseMatrix))
            throw new UnsupportedOperationException("X must be a DenseMatrix");

        checkSolve(B, X);

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        // Allocate factorization matrix. The factorization matrix will be
        // large enough to accomodate any pivots
        BandMatrix Af = new BandMatrix(this, kd, kd + kd);
        int[] ipiv = new int[numRows];

        intW info = new intW(0);
        LAPACK.getInstance().dgbsv(numRows, kd, kd, X.numColumns(),
                Af.getData(), Matrices.ld(2 * kd + kd + 1), ipiv, Xd,
                Matrices.ld(numRows), info);

        if (info.val > 0)
            throw new MatrixSingularException();
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return X;
    }

    @Override
    public Vector solve(Vector b, Vector x) {
        DenseMatrix B = new DenseMatrix(b, false), X = new DenseMatrix(x, false);
        solve(B, X);
        return x;
    }

    @Override
    public Matrix transSolve(Matrix B, Matrix X) {
        return solve(B, X);
    }

    @Override
    public Vector transSolve(Vector b, Vector x) {
        return solve(b, x);
    }

    Matrix SPDsolve(Matrix B, Matrix X) {
        if (!(X instanceof DenseMatrix))
            throw new UnsupportedOperationException("X must be a DenseMatrix");

        checkSolve(B, X);

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        intW info = new intW(0);
        LAPACK.getInstance().dpbsv(uplo.netlib(), numRows, kd, X.numColumns(),
                data.clone(), Matrices.ld(kd + 1), Xd, Matrices.ld(numRows), info);

        if (info.val > 0)
            throw new MatrixNotSPDException();
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return X;
    }

    @Override
    public Matrix transpose() {
        return this;
    }

}
