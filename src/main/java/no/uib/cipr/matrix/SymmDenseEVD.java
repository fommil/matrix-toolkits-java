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
import org.netlib.util.intW;


/**
 * Computes eigenvalues of symmetrical, dense matrices
 */
public class SymmDenseEVD extends SymmEVD {

    /**
     * Double work array
     */
    private final double[] work;

    /**
     * Integer work array
     */
    private final int[] iwork;

    /**
     * Upper or lower part stored
     */
    private final UpLo uplo;

    /**
     * Range of eigenvalues to compute
     */
    private final JobEigRange range;

    /**
     * Eigenvector supports
     */
    private final int[] isuppz;

    /**
     * Tolerance criteria
     */
    private final double abstol;

    /**
     * Sets up an eigenvalue decomposition for symmetrical, dense matrices.
     * Computes all eigenvalues and eigenvectors, and uses a low default
     * tolerance criteria
     * 
     * @param n
     *            Size of the matrix
     * @param upper
     *            True if the upper part of the matrix is stored, and false if
     *            the lower part of the matrix is stored instead
     */
    public SymmDenseEVD(int n, boolean upper) {
        this(n, upper, true, LAPACK.getInstance().dlamch("Safe minimum"));
    }

    /**
     * Sets up an eigenvalue decomposition for symmetrical, dense matrices.
     * Computes all eigenvalues and eigenvectors
     * 
     * @param n
     *            Size of the matrix
     * @param upper
     *            True if the upper part of the matrix is stored, and false if
     *            the lower part of the matrix is stored instead
     * @param abstol
     *            Absolute tolerance criteria
     */
    public SymmDenseEVD(int n, boolean upper, double abstol) {
        this(n, upper, true, abstol);
    }

    /**
     * Sets up an eigenvalue decomposition for symmetrical, dense matrices. Uses
     * a low default tolerance criteria
     * 
     * @param n
     *            Size of the matrix
     * @param upper
     *            True if the upper part of the matrix is stored, and false if
     *            the lower part of the matrix is stored instead
     * @param vectors
     *            True to compute the eigenvectors, false for just the
     *            eigenvalues
     */
    public SymmDenseEVD(int n, boolean upper, boolean vectors) {
        this(n, upper, vectors, LAPACK.getInstance().dlamch("Safe minimum"));
    }

    /**
     * Sets up an eigenvalue decomposition for symmetrical, dense matrices
     * 
     * @param n
     *            Size of the matrix
     * @param upper
     *            True if the upper part of the matrix is stored, and false if
     *            the lower part of the matrix is stored instead
     * @param vectors
     *            True to compute the eigenvectors, false for just the
     *            eigenvalues
     * @param abstol
     *            Absolute tolerance criteria
     */
    public SymmDenseEVD(int n, boolean upper, boolean vectors, double abstol) {
        super(n, vectors);
        this.abstol = abstol;

        uplo = upper ? UpLo.Upper : UpLo.Lower;
        range = JobEigRange.All;
        isuppz = new int[2 * Math.max(1, n)];

        // Find the needed workspace
        double[] worksize = new double[1];
        int[] iworksize = new int[1];
        intW info = new intW(0);
        LAPACK.getInstance().dsyevr(job.netlib(), range.netlib(), uplo.netlib(), n,
        	new double[0], Matrices.ld(n), 0, 0, 0, 0, abstol, new intW(1), new double[0], new double[0],
        	Matrices.ld(n), isuppz, worksize, -1, iworksize, -1, info);

        // Allocate workspace
        int lwork = 0, liwork = 0;
        if (info.val != 0) {
            lwork = 26 * n;
            liwork = 10 * n;
        } else {
            lwork = (int) worksize[0];
            liwork = iworksize[0];
        }

        lwork = Math.max(1, lwork);
        liwork = Math.max(1, liwork);
        work = new double[lwork];
        iwork = new int[liwork];
    }

    /**
     * Convenience method for computing the full eigenvalue decomposition of the
     * given matrix
     * 
     * @param A
     *            Matrix to factorize. Upper part extracted, and the matrix is
     *            not modified
     * @return Newly allocated decomposition
     * @throws NotConvergedException
     */
    public static SymmDenseEVD factorize(Matrix A) throws NotConvergedException {
        return new SymmDenseEVD(A.numRows(), true)
                .factor(new UpperSymmDenseMatrix(A));
    }

    /**
     * Computes the eigenvalue decomposition of the given matrix
     * 
     * @param A
     *            Matrix to factorize. Overwritten on return
     * @return The current eigenvalue decomposition
     * @throws NotConvergedException
     */
    public SymmDenseEVD factor(LowerSymmDenseMatrix A)
            throws NotConvergedException {
        if (uplo != UpLo.Lower)
            throw new IllegalArgumentException(
                    "Eigenvalue computer configured for lower-symmetrical matrices");

        return factor(A, A.getData());
    }

    /**
     * Computes the eigenvalue decomposition of the given matrix
     * 
     * @param A
     *            Matrix to factorize. Overwritten on return
     * @return The current eigenvalue decomposition
     * @throws NotConvergedException
     */
    public SymmDenseEVD factor(UpperSymmDenseMatrix A)
            throws NotConvergedException {
        if (uplo != UpLo.Upper)
            throw new IllegalArgumentException(
                    "Eigenvalue computer configured for upper-symmetrical matrices");

        return factor(A, A.getData());
    }

    private SymmDenseEVD factor(Matrix A, double[] data)
            throws NotConvergedException {
        if (A.numRows() != n)
            throw new IllegalArgumentException("A.numRows() != n");

        intW info = new intW(0);
        LAPACK.getInstance().dsyevr(job.netlib(), range.netlib(), uplo.netlib(), n, data,
        	Matrices.ld(n), 0, 0, 0, 0, abstol, new intW(1), w,
                job == JobEig.All ? Z.getData() : new double[0], Matrices.ld(n), isuppz, work,
                work.length, iwork, iwork.length, info);

        if (info.val > 0)
            throw new NotConvergedException(
                    NotConvergedException.Reason.Iterations);
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return this;
    }

}
