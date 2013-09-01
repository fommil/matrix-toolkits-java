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
 * Computes eigenvalue decompositions of general matrices
 */
public class EVD {

    /**
     * Double work array
     */
    private final double[] work;

    /**
     * Size of the matrix
     */
    private final int n;

    /**
     * Job to do on the left and right eigenvectors
     */
    private final JobEig jobLeft, jobRight;

    /**
     * Contains the real and imaginary parts of the eigenvalues
     */
    private final double[] Wr, Wi;

    /**
     * Contains the left and the right eigenvectors
     */
    private final DenseMatrix Vl, Vr;

    /**
     * Creates an empty eigenvalue decomposition which will compute all the
     * eigenvalues and eigenvectors (left and right)
     * 
     * @param n
     *            Size of the matrix
     */
    public EVD(int n) {
        this(n, true, true);
    }

    /**
     * Creates an empty eigenvalue decomposition
     * 
     * @param n
     *            Size of the matrix
     * @param left
     *            Whether to compute the left eigenvectors or not
     * @param right
     *            Whether to compute the right eigenvectors or not
     */
    public EVD(int n, boolean left, boolean right) {
        this.n = n;
        this.jobLeft = left ? JobEig.All : JobEig.Eigenvalues;
        this.jobRight = right ? JobEig.All : JobEig.Eigenvalues;

        // Allocate space for the decomposition
        Wr = new double[n];
        Wi = new double[n];

        if (left)
            Vl = new DenseMatrix(n, n);
        else
            Vl = null;

        if (right)
            Vr = new DenseMatrix(n, n);
        else
            Vr = null;

        // Find the needed workspace
        double[] worksize = new double[1];
        intW info = new intW(0);
        LAPACK.getInstance().dgeev(jobLeft.netlib(), jobRight.netlib(), n, new double[0],
                Matrices.ld(n), new double[0], new double[0], new double[0], Matrices.ld(n),
                new double[0], Matrices.ld(n), worksize, -1, info);

        // Allocate workspace
        int lwork = 0;
        if (info.val != 0) {
            if (jobLeft == JobEig.All || jobRight == JobEig.All)
                lwork = 4 * n;
            else
                lwork = 3 * n;
        } else
            lwork = (int) worksize[0];

        lwork = Math.max(1, lwork);
        work = new double[lwork];
    }

    /**
     * Convenience method for computing the complete eigenvalue decomposition of
     * the given matrix
     * 
     * @param A
     *            Matrix to factorize. Not modified
     * @return Newly allocated decomposition
     * @throws NotConvergedException
     */
    public static EVD factorize(Matrix A) throws NotConvergedException {
        return new EVD(A.numRows()).factor(new DenseMatrix(A));
    }

    /**
     * Computes the eigenvalue decomposition of the given matrix
     * 
     * @param A
     *            Matrix to factorize. Overwritten on return
     * @return The current decomposition
     * @throws NotConvergedException
     */
    public EVD factor(DenseMatrix A) throws NotConvergedException {
        if (!A.isSquare())
            throw new IllegalArgumentException("!A.isSquare()");
        else if (A.numRows() != n)
            throw new IllegalArgumentException("A.numRows() != n");

        intW info = new intW(0);
        LAPACK.getInstance().dgeev(jobLeft.netlib(), jobRight.netlib(), n, A.getData(),
                Matrices.ld(n), Wr, Wi, jobLeft == JobEig.All ? Vl.getData() : new double[0],
                Matrices.ld(n), jobRight == JobEig.All ? Vr.getData() : new double[0], Matrices.ld(n),
                work, work.length, info);

        if (info.val > 0)
            throw new NotConvergedException(
                    NotConvergedException.Reason.Iterations);
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return this;
    }

    /**
     * Gets the left eigenvectors, if available
     */
    public DenseMatrix getLeftEigenvectors() {
        return Vl;
    }

    /**
     * Gets the right eigenvectors, if available
     */
    public DenseMatrix getRightEigenvectors() {
        return Vr;
    }

    /**
     * Gets the real part of the eigenvalues
     */
    public double[] getRealEigenvalues() {
        return Wr;
    }

    /**
     * Gets the imaginary part of the eigenvalues
     */
    public double[] getImaginaryEigenvalues() {
        return Wi;
    }

    /**
     * True if the left eigenvectors have been computed
     */
    public boolean hasLeftEigenvectors() {
        return Vl != null;
    }

    /**
     * True if the right eigenvectors have been computed
     */
    public boolean hasRightEigenvectors() {
        return Vr != null;
    }

}
