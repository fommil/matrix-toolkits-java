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
 * Partial implementation of a symmetrical, dense matrix
 */
abstract class AbstractSymmDenseMatrix extends AbstractDenseMatrix {

    /**
     * Upper or lower part stored?
     */
    private UpLo uplo;

    /**
     * Constructor for AbstractSymmDenseMatrix
     */
    AbstractSymmDenseMatrix(int n, UpLo uplo) {
        super(n, n);
        this.uplo = uplo;
    }

    /**
     * Constructor for AbstractSymmDenseMatrix
     */
    AbstractSymmDenseMatrix(Matrix A, UpLo uplo) {
        this(A, true, uplo);
    }

    /**
     * Constructor for AbstractSymmDenseMatrix
     */
    AbstractSymmDenseMatrix(Matrix A, boolean deep, UpLo uplo) {
        super(A, deep);
        if (!isSquare())
            throw new IllegalArgumentException(
                    "Symmetric matrix must be square");
        this.uplo = uplo;
    }

    @Override
    public Matrix multAdd(double alpha, Matrix B, Matrix C) {
        if (!(B instanceof DenseMatrix) || !(C instanceof DenseMatrix))
            return super.multAdd(alpha, B, C);

        checkMultAdd(B, C);

        double[] Bd = ((DenseMatrix) B).getData(), Cd = ((DenseMatrix) C)
                .getData();

        BLAS.getInstance().dsymm(Side.Left.netlib(), uplo.netlib(), C.numRows(), C.numColumns(),
                alpha, data, Math.max(1, C.numRows()), Bd,
                Math.max(1, C.numRows()), 1, Cd, Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public Matrix transAmultAdd(double alpha, Matrix B, Matrix C) {
        return multAdd(alpha, B, C);
    }

    @Override
    public Matrix rank1(double alpha, Vector x, Vector y) {
        if (x != y)
            throw new IllegalArgumentException("x != y");
        if (!(x instanceof DenseVector))
            return super.rank1(alpha, x, y);

        checkRank1(x, y);

        double[] xd = ((DenseVector) x).getData();

        BLAS.getInstance().dsyr(uplo.netlib(), numRows, alpha, xd, 1, data,
                Math.max(1, numRows));

        return this;
    }

    @Override
    public Matrix rank2(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.rank2(alpha, x, y);

        checkRank2(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dsyr2(uplo.netlib(), numRows, alpha, xd, 1, yd, 1, data,
                Math.max(1, numRows));

        return this;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.multAdd(alpha, x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dsymv(uplo.netlib(), numRows, alpha, data, Math.max(1, numRows),
                xd, 1, 1, yd, 1);

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        return multAdd(alpha, x, y);
    }

    @Override
    public Matrix rank1(double alpha, Matrix C) {
        if (!(C instanceof DenseMatrix))
            return super.rank1(alpha, C);

        checkRank1(C);

        double[] Cd = ((DenseMatrix) C).getData();

        BLAS.getInstance().dsyrk(uplo.netlib(), Transpose.NoTranspose.netlib(), numRows,
                C.numColumns(), alpha, Cd, Math.max(1, numRows), 1, data,
                Math.max(1, numRows));

        return this;
    }

    @Override
    public Matrix transRank1(double alpha, Matrix C) {
        if (!(C instanceof DenseMatrix))
            return super.transRank1(alpha, C);

        checkTransRank1(C);

        double[] Cd = ((DenseMatrix) C).getData();

        BLAS.getInstance().dsyrk(uplo.netlib(), Transpose.Transpose.netlib(), numRows, numRows,
                alpha, Cd, Math.max(1, numRows), 1, data, Math.max(1, numRows));

        return this;
    }

    @Override
    public Matrix rank2(double alpha, Matrix B, Matrix C) {
        if (!(B instanceof DenseMatrix) || !(C instanceof DenseMatrix))
            return super.rank2(alpha, B, C);

        checkRank2(B, C);

        double[] Bd = ((DenseMatrix) B).getData(), Cd = ((DenseMatrix) C)
                .getData();

        BLAS.getInstance().dsyr2k(uplo.netlib(), Transpose.NoTranspose.netlib(), numRows,
                B.numColumns(), alpha, Bd, Math.max(1, numRows), Cd,
                Math.max(1, numRows), 1, data, Math.max(1, numRows));

        return this;
    }

    @Override
    public Matrix transRank2(double alpha, Matrix B, Matrix C) {
        if (!(B instanceof DenseMatrix) || !(C instanceof DenseMatrix))
            return super.transRank2(alpha, B, C);

        checkTransRank2(B, C);

        double[] Bd = ((DenseMatrix) B).getData(), Cd = ((DenseMatrix) C)
                .getData();

        BLAS.getInstance().dsyr2k(uplo.netlib(), Transpose.Transpose.netlib(), numRows, B.numRows(),
                alpha, Bd, Math.max(1, B.numRows()), Cd,
                Math.max(1, B.numRows()), 1, data, Math.max(1, numRows));

        return this;
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        if (!(X instanceof DenseMatrix))
            throw new UnsupportedOperationException("X must be a DenseMatrix");

        checkSolve(B, X);

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        // Allocate factorization matrix
        double[] newData = data.clone();
        int[] ipiv = new int[numRows];

        // Query optimal workspace
        double[] work = new double[1];
        intW info = new intW(0);
        LAPACK.getInstance().dsysv(uplo.netlib(), numRows, X.numColumns(),
                newData, Matrices.ld(numRows), ipiv, Xd, Matrices.ld(numRows),
                work, -1, info);

        // Allocate workspace
        int lwork = -1;
        if (info.val != 0)
            lwork = 1;
        else
            lwork = Math.max((int) work[0], 1);
        work = new double[lwork];

        // Solve
        info.val = 0;
        LAPACK.getInstance().dsysv(uplo.netlib(), numRows, X.numColumns(), newData,
                Matrices.ld(numRows), ipiv, Xd, Matrices.ld(numRows), work, lwork, info);

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
        LAPACK.getInstance().dposv(uplo.netlib(), numRows, X.numColumns(),
                data.clone(), Matrices.ld(numRows), Xd, Matrices.ld(numRows), info);

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
