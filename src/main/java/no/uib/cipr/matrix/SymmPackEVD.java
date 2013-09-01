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
 * Computes eigenvalues of symmetrical, packed matrices
 */
public class SymmPackEVD extends SymmEVD {

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
     * Sets up an eigenvalue decomposition for symmetrical, packed matrices.
     * Computes all eigenvalues and eigenvectors
     * 
     * @param n
     *            Size of the matrix
     * @param upper
     *            True if the upper part of the matrix is stored, and false if
     *            the lower part of the matrix is stored instead
     */
    public SymmPackEVD(int n, boolean upper) {
        this(n, upper, true);
    }

    /**
     * Sets up an eigenvalue decomposition for symmetrical, packed matrices
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
    public SymmPackEVD(int n, boolean upper, boolean vectors) {
        super(n, vectors);

        uplo = upper ? UpLo.Upper : UpLo.Lower;

        // Find the needed workspace
        double[] worksize = new double[1];
        int[] iworksize = new int[1];
        intW info = new intW(0);
        LAPACK.getInstance().dspevd(job.netlib(), uplo.netlib(), n, new double[0],
                new double[0], new double[0], Matrices.ld(n), worksize, -1, iworksize, -1, info);

        // Allocate workspace
        int lwork = 0, liwork = 0;
        if (info.val != 0) {
            if (job == JobEig.All) {
                lwork = 1 + 6 * n + n * n;
                liwork = 3 + 5 * n;
            } else {
                lwork = 2 * n;
                liwork = 1;
            }
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
    public static SymmPackEVD factorize(Matrix A) throws NotConvergedException {
        return new SymmPackEVD(A.numRows(), true)
                .factor(new UpperSymmPackMatrix(A));
    }

    /**
     * Computes the eigenvalue decomposition of the given matrix
     * 
     * @param A
     *            Matrix to factorize. Overwritten on return
     * @return The current eigenvalue decomposition
     * @throws NotConvergedException
     */
    public SymmPackEVD factor(LowerSymmPackMatrix A)
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
    public SymmPackEVD factor(UpperSymmPackMatrix A)
            throws NotConvergedException {
        if (uplo != UpLo.Upper)
            throw new IllegalArgumentException(
                    "Eigenvalue computer configured for upper-symmetrical matrices");

        return factor(A, A.getData());
    }

    private SymmPackEVD factor(Matrix A, double[] data)
            throws NotConvergedException {
        if (A.numRows() != n)
            throw new IllegalArgumentException("A.numRows() != n");

        intW info = new intW(0);
        LAPACK.getInstance().dspevd(job.netlib(), uplo.netlib(), n, data, w,
                job == JobEig.All ? Z.getData() : new double[0], Matrices.ld(n), work,
                work.length, iwork, iwork.length, info);

        if (info.val > 0)
            throw new NotConvergedException(
                    NotConvergedException.Reason.Iterations);
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return this;
    }

}
