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

import no.uib.cipr.matrix.Matrix.Norm;

import com.github.fommil.netlib.LAPACK;
import org.netlib.util.doubleW;
import org.netlib.util.intW;

/**
 * Dense Cholesky decomposition
 */
public class DenseCholesky {

    /**
     * Matrix dimension
     */
    private final int n;

    /**
     * Cholesky decomposition of a lower matrix
     */
    private LowerTriangDenseMatrix Cl;

    /**
     * Cholesky decomposition of an upper matrix
     */
    private UpperTriangDenseMatrix Cu;

    /**
     * If the matrix is SPD or not
     */
    private boolean notspd;

    /**
     * True for upper part, else false
     */
    private final boolean upper;

    /**
     * Constructor for DenseCholesky
     * 
     * @param n
     *            Matrix size
     * @param upper
     *            True for decomposing an upper symmetrical matrix, false for a
     *            lower symmetrical matrix
     */
    public DenseCholesky(int n, boolean upper) {
        this.n = n;
        this.upper = upper;

        if (upper)
            Cu = new UpperTriangDenseMatrix(n);
        else
            Cl = new LowerTriangDenseMatrix(n);
    }

    /**
     * Calculates a Cholesky decomposition
     * 
     * @param A
     *            Matrix to decompose. Not modified
     * @return The current decomposition
     */
    public static DenseCholesky factorize(Matrix A) {
        return new DenseCholesky(A.numRows(), true)
                .factor(new UpperSPDDenseMatrix(A));
    }

    /**
     * Calculates a Cholesky decomposition
     * 
     * @param A
     *            Matrix to decompose. Overwritten on return
     * @return The current decomposition
     */
    public DenseCholesky factor(LowerSPDDenseMatrix A) {
        if (upper)
            throw new IllegalArgumentException(
                    "Cholesky decomposition constructed for upper matrices");

        return decompose(A);
    }

    /**
     * Calculates a Cholesky decomposition
     * 
     * @param A
     *            Matrix to decompose. Overwritten on return
     * @return The current decomposition
     */
    public DenseCholesky factor(UpperSPDDenseMatrix A) {
        if (!upper)
            throw new IllegalArgumentException(
                    "Cholesky decomposition constructed for lower matrices");

        return decompose(A);
    }

    private DenseCholesky decompose(AbstractDenseMatrix A) {
        if (n != A.numRows())
            throw new IllegalArgumentException("n != A.numRows()");

        notspd = false;

        intW info = new intW(0);
        if (upper)
            LAPACK.getInstance().dpotrf(UpLo.Upper.netlib(), A.numRows(),
                    A.getData(), Matrices.ld(A.numRows()), info);
        else
            LAPACK.getInstance().dpotrf(UpLo.Lower.netlib(), A.numRows(),
                    A.getData(), Matrices.ld(A.numRows()), info);

        if (info.val > 0)
            notspd = true;
        else if (info.val < 0)
            throw new IllegalArgumentException();

        if (upper)
            Cu.set(A);
        else
            Cl.set(A);

        return this;
    }

    /**
     * Returns true if the matrix decomposed is symmetrical, positive definite
     */
    public boolean isSPD() {
        return !notspd;
    }

    /**
     * Returns the decomposition matrix. Only valid for decomposition of a lower
     * SPD matrix
     */
    public LowerTriangDenseMatrix getL() {
        if (!upper)
            return Cl;
        else
            throw new UnsupportedOperationException();
    }

    /**
     * Returns the decomposition matrix. Only valid for decomposition of a upper
     * SPD matrix
     */
    public UpperTriangDenseMatrix getU() {
        if (upper)
            return Cu;
        else
            throw new UnsupportedOperationException();
    }

    /**
     * Solves for <code>B</code>, overwriting it on return
     */
    public DenseMatrix solve(DenseMatrix B) throws MatrixNotSPDException {
        if (notspd)
            throw new MatrixNotSPDException();
        if (n != B.numRows())
            throw new IllegalArgumentException("n != B.numRows()");

        intW info = new intW(0);
        if (upper)
            LAPACK.getInstance().dpotrs(UpLo.Upper.netlib(), Cu.numRows(),
                    B.numColumns(), Cu.getData(), Matrices.ld(Cu.numRows()),
                    B.getData(), Matrices.ld(Cu.numRows()), info);
        else
            LAPACK.getInstance().dpotrs(UpLo.Lower.netlib(), Cl.numRows(),
                    B.numColumns(), Cl.getData(), Matrices.ld(Cl.numRows()),
                    B.getData(), Matrices.ld(Cl.numRows()), info);

        if (info.val < 0)
            throw new IllegalArgumentException();

        return B;
    }

    /**
     * Computes the reciprocal condition number
     * 
     * @param A
     *            The matrix this is a decomposition of
     * @return The reciprocal condition number. Values close to unity indicate a
     *         well-conditioned system, while numbers close to zero do not.
     */
    public double rcond(Matrix A) {
        if (n != A.numRows())
            throw new IllegalArgumentException("n != A.numRows()");
        if (!A.isSquare())
            throw new IllegalArgumentException("!A.isSquare()");

        double anorm = A.norm(Norm.One);

        double[] work = new double[3 * n];
        int[] iwork = new int[n];

        intW info = new intW(0);
        doubleW rcond = new doubleW(0);
        if (upper)
            LAPACK.getInstance().dpocon(UpLo.Upper.netlib(), n, Cu.getData(),
            	Matrices.ld(n), anorm, rcond, work, iwork, info);
        else
            LAPACK.getInstance().dpocon(UpLo.Lower.netlib(), n, Cl.getData(),
            	Matrices.ld(n), anorm, rcond, work, iwork, info);

        if (info.val < 0)
            throw new IllegalArgumentException();

        return rcond.val;
    }

}
