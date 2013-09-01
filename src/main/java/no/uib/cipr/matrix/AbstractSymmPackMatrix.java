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
 * Partial implementation of a symmetrical, packed matrix
 */
abstract class AbstractSymmPackMatrix extends AbstractPackMatrix {

    /**
     * Which part of the matrix which is stored
     */
    private UpLo uplo;

    /**
     * Constructor for AbstractSymmPackMatrix
     */
    AbstractSymmPackMatrix(int n, UpLo uplo) {
        super(n);
        this.uplo = uplo;
    }

    /**
     * Constructor for AbstractSymmPackMatrix
     */
    AbstractSymmPackMatrix(Matrix A, UpLo uplo) {
        this(A, true, uplo);
    }

    /**
     * Constructor for AbstractSymmPackMatrix
     */
    AbstractSymmPackMatrix(Matrix A, boolean deep, UpLo uplo) {
        super(A, deep);
        this.uplo = uplo;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.multAdd(alpha, x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dspmv(uplo.netlib(), numRows, alpha, data, xd, 1, 1, yd, 1);

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        return multAdd(alpha, x, y);
    }

    @Override
    public Matrix rank1(double alpha, Vector x, Vector y) {
        if (x != y)
            throw new IllegalArgumentException("x != y");
        if (!(x instanceof DenseVector))
            return super.rank1(alpha, x, y);

        checkRank1(x, y);

        double[] xd = ((DenseVector) x).getData();

        BLAS.getInstance().dspr(uplo.netlib(), numRows, alpha, xd, 1, data);

        return this;
    }

    @Override
    public Matrix rank2(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.rank2(alpha, x, y);

        checkRank2(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dspr2(uplo.netlib(), numRows, alpha, xd, 1, yd, 1, data);

        return this;
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        if (!(X instanceof DenseMatrix))
            throw new UnsupportedOperationException("X must be a DenseMatrix");

        checkSolve(B, X);

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        int[] ipiv = new int[numRows];

        intW info = new intW(0);
        LAPACK.getInstance().dspsv(uplo.netlib(), numRows, X.numColumns(),
                data.clone(), ipiv, Xd, Matrices.ld(numRows), info);

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
        LAPACK.getInstance().dppsv(uplo.netlib(), numRows, X.numColumns(),
                data.clone(), Xd, Matrices.ld(numRows), info);

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
