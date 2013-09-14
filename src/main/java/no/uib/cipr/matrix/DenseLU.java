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

import com.github.fommil.netlib.LAPACK;
import no.uib.cipr.matrix.Matrix.Norm;
import org.netlib.util.doubleW;
import org.netlib.util.intW;

/**
 * Dense Partial Pivot LU decomposition: {@code A = P * L * U}.
 */
public class DenseLU {

    /**
     * Holds the LU factors
     */
    private DenseMatrix LU;

    /**
     * Row pivotations
     */
    private int[] piv;

    /**
     * True if the matrix was singular
     */
    private boolean singular;

    /**
     * Constructor for DenseLU
     * 
     * @param m
     *            Number of rows
     * @param n
     *            Number of columns
     */
    public DenseLU(int m, int n) {
        LU = new DenseMatrix(m, n);
        piv = new int[Math.min(m, n)];
    }

    /**
     * Creates an LU decomposition of the given matrix
     * 
     * @param A
     *            Matrix to decompose. Not modified
     * @return The current decomposition
     */
    public static DenseLU factorize(Matrix A) {
        return new DenseLU(A.numRows(), A.numColumns()).factor(new DenseMatrix(A));
    }

    /**
     * Creates an LU decomposition of the given matrix
     * 
     * @param A
     *            Matrix to decompose. Overwritten with the decomposition
     * @return The current decomposition
     */
    public DenseLU factor(DenseMatrix A) {
        singular = false;

        intW info = new intW(0);
        LAPACK.getInstance().dgetrf(A.numRows(), A.numColumns(),
                A.getData(), Matrices.ld(A.numRows()), piv, info);

        if (info.val > 0)
            singular = true;
        else if (info.val < 0)
            throw new IllegalArgumentException();

        LU.set(A);

        return this;
    }

    /**
     * Returns the permutation matrix.
     */
    public PermutationMatrix getP() {
      PermutationMatrix perm = PermutationMatrix.fromPartialPivots(piv);
      perm.transpose();
      return perm;
    }

    /**
     * Returns the lower triangular factor
     */
    public UnitLowerTriangDenseMatrix getL() {
        return new UnitLowerTriangDenseMatrix(getLU(), false);
    }

    /**
     * Returns the upper triangular factor
     */
    public UpperTriangDenseMatrix getU() {
        return new UpperTriangDenseMatrix(getLU(), false);
    }

    /**
     * Returns the decomposition matrix
     */
    protected DenseMatrix getLU() {
        return LU;
    }

    /**
     * Computes the reciprocal condition number, using either the infinity norm
     * of the 1 norm.
     * 
     * @param A
     *            The matrix this is a decomposition of
     * @param norm
     *            Either <code>Norm.One</code> or <code>Norm.Infinity</code>
     * @return The reciprocal condition number. Values close to unity indicate a
     *         well-conditioned system, while numbers close to zero do not.
     */
    public double rcond(Matrix A, Norm norm) {
        if (norm != Norm.One && norm != Norm.Infinity)
            throw new IllegalArgumentException(
                    "Only the 1 or the Infinity norms are supported");

        double anorm = A.norm(norm);

        int n = A.numRows();

        intW info = new intW(0);
        doubleW rcond = new doubleW(0);
        LAPACK.getInstance().dgecon(norm.netlib(), n, LU.getData(), Matrices.ld(n), anorm,
                rcond, new double[4 * n], new int[n], info);

        if (info.val < 0)
            throw new IllegalArgumentException();

        return rcond.val;
    }

    /**
     * Returns the row pivots
     */
    public int[] getPivots() {
        return piv;
    }

    /**
     * Checks for singularity
     */
    public boolean isSingular() {
        return singular;
    }

    /**
     * Computes <code>A\B</code>, overwriting <code>B</code>
     */
    public DenseMatrix solve(DenseMatrix B) throws MatrixSingularException {
        return solve(B, Transpose.NoTranspose);
    }

    /**
     * Computes <code>A<sup>T</sup>\B</code>, overwriting <code>B</code>
     */
    public DenseMatrix transSolve(DenseMatrix B) throws MatrixSingularException {
        return solve(B, Transpose.Transpose);
    }

    private DenseMatrix solve(DenseMatrix B, Transpose trans)
            throws MatrixSingularException {
        if (singular)
            throw new MatrixSingularException();
        if (B.numRows() != LU.numRows())
            throw new IllegalArgumentException("B.numRows() != LU.numRows()");

        intW info = new intW(0);
        LAPACK.getInstance().dgetrs(trans.netlib(), LU.numRows(),
                B.numColumns(), LU.getData(), Matrices.ld(LU.numRows()),
                piv, B.getData(), Matrices.ld(LU.numRows()), info);

        if (info.val < 0)
            throw new IllegalArgumentException();

        return B;
    }

}
