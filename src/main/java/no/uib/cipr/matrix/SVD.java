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
 * Computes singular value decompositions
 */
public class SVD {

    /**
     * Work array
     */
    private final double[] work;

    /**
     * Work array
     */
    private final int[] iwork;

    /**
     * Matrix dimension
     */
    private final int m, n;

    /**
     * Job to do. <br>
     * if job2 = null then GDESDD is used else GDESVD is used.
     */
    private final JobSVD job1, job2;

    /**
     * The singular values
     */
    private final double[] S;

    /**
     * Singular vectors
     */
    private DenseMatrix U, Vt;

    /**
     * Leading Dimensions of A, U, and Vt.
     */
    private final int lda, ldu, ldvt;

    /**
     * Creates an empty SVD which will compute all singular values and vectors
     * 
     * @param m
     *            Number of rows
     * @param n
     *            Number of columns
     */
    public SVD(int m, int n) {
        this(m, n, true);
    }

    /**
     * Creates an empty SVD
     * 
     * @param m
     *            Number of rows
     * @param n
     *            Number of columns
     * @param vectors
     *            True to compute the singular vectors, false for just the
     *            singular values
     */
    public SVD(int m, int n, boolean vectors) {
        this(m, n, vectors ? JobSVD.All : JobSVD.None);
    }

    /**
     * Creates and empty SVD which will apply the {@link LAPACK#dgesvd}
     * algorithm. Compared to {@link LAPACK#dgesdd}, this is generally slower,
     * uses less memory and is more numerically stable on some implementations.<br>
     * {@link JobSVD} arguments <code>jobU</code> and <code>jobVT</code>
     * independently control the left and right singular vector computation.<br>
     * 
     * jobU and jobVT cannot both be Overwrite.
     * 
     * @param m
     *            Number of rows
     * @param n
     *            Number of Columns
     * @param jobU
     *            Job for left singular vectors, (All|None|Partial|Overwrite)
     * @param jobVT
     *            Job for right singular vectors, (All|None|Partial|Overwrite)
     */
    public SVD(int m, int n, JobSVD jobU, JobSVD jobVT) {
        this.m = m;
        this.n = n;
        int min = Math.min(m, n);
        job1 = jobU;
        job2 = jobVT;

        lda = Matrices.ld(m);

        if (jobU == JobSVD.None) { // N
            ldu = 1;
            U = new DenseMatrix(0, 0);
        } else if (jobU == JobSVD.Part) { // S
            ldu = m;
            U = new DenseMatrix(ldu, min);
        } else if (jobU == JobSVD.Overwrite) { // O
            ldu = 1;
            U = new DenseMatrix(0, 0);
        } else { // All //A
            ldu = m;
            U = new DenseMatrix(ldu, m);
        }

        if (jobVT == JobSVD.None) { // N
            ldvt = 1;
            Vt = new DenseMatrix(0, 0);
        } else if (jobVT == JobSVD.Part) { // S
            ldvt = min;
            Vt = new DenseMatrix(ldvt, n);
        } else if (jobVT == JobSVD.Overwrite) { // O
            ldvt = 1;
            Vt = new DenseMatrix(0, 0);
        } else { // All
            ldvt = n;
            Vt = new DenseMatrix(ldvt, n);
        }

        int lwork = queryOptimalDgesvdWorksize(m, n, lda, ldu, ldvt, jobU,
                jobVT);

        // Allocate space for the decomposition
        S = new double[min];
        iwork = null;
        work = new double[Math.max(1, lwork)];
    }

    /**
     * Creates and empty SVD which will apply the {@link LAPACK#dgesdd}
     * algorithm.<br>
     * Compared to the {@link LAPACK#dgesvd} algorithm, this is typically 5-10x
     * faster for large matrices, potentially uses more memory, and some
     * implementations are less numerically stable.<br>
     * <br>
     * When {@link JobSVD} argument <code>jobz</code> is Overwrite:<br>
     * if m >= n then<br>
     * &emsp;first m cols of U and written to A<br>
     * else<br>
     * &emsp;first n rows of Vt are written to A<br>
     * 
     * @param m
     * @param n
     * @param jobz
     *            Job for both singular vectors, (All|None|Partial|Overwrite)
     */
    public SVD(int m, int n, JobSVD jobz) {
        this.m = m;
        this.n = n;
        int min = Math.min(m, n);
        job1 = jobz;
        job2 = null;

        lda = Matrices.ld(m);

        // Find workspace requirements
        if (jobz == JobSVD.None) { // N
            ldu = 1;
            U = new DenseMatrix(0, 0); // ldu
            ldvt = 1;
            Vt = new DenseMatrix(0, 0); // ldvt
        } else if (jobz == JobSVD.Part) { // S
            ldu = Matrices.ld(m);
            U = new DenseMatrix(ldu, min);
            ldvt = Matrices.ld(min);
            Vt = new DenseMatrix(ldvt, n);
        } else if (jobz == JobSVD.Overwrite) { // O
            if (m >= n) {
                ldvt = Matrices.ld(n);
                Vt = new DenseMatrix(ldvt, n);
                ldu = 1;
                U = new DenseMatrix(0, 0); // ldu
            } else {
                ldvt = 1;
                Vt = new DenseMatrix(0, 0); // ldvt
                ldu = Matrices.ld(m);
                U = new DenseMatrix(ldu, m);
            }
        } else { // All //A
            ldu = Matrices.ld(m);
            U = new DenseMatrix(ldu, m);
            ldvt = Matrices.ld(n);
            Vt = new DenseMatrix(ldvt, n);
        }

        int lwork = queryOptimalDgesddWorksize(m, n, lda, ldu, ldvt, jobz);

        // Allocate space for the decomposition
        S = new double[min];
        iwork = new int[8 * Math.min(m, n)];
        work = new double[lwork];
    }

    // Query optimal workspace
    private static int queryOptimalDgesvdWorksize(int m, int n, int lda,
            int ldu, int ldvt, JobSVD jobU, JobSVD jobVT) {
        int min = Math.min(m, n);
        int max = Math.max(m, n);
        double[] worksize = new double[1];
        intW info = new intW(0);
        LAPACK.getInstance().dgesvd(jobU.netlib(), jobVT.netlib(), m, n,
                new double[0], lda, // A
                new double[0], // S
                new double[0], ldu, // U
                new double[0], ldvt, // Vt
                worksize, -1, info);

        int lwork = Math.max(3 * min + max, 5 * min);
        if (info.val == 0) {
            lwork = (int) worksize[0];
        }

        return Math.max(lwork, 1);
    }

    // Query optimal workspace
    private static int queryOptimalDgesddWorksize(int m, int n, int lda,
            int ldu, int ldvt, JobSVD jobz) {
        int min = Math.min(m, n);
        int max = Math.max(m, n);

        // Query optimal workspace
        double[] worksize = new double[1];
        intW info = new intW(0);
        LAPACK.getInstance().dgesdd(jobz.netlib(), m, n, new double[0], lda, // A
                new double[0], // S
                new double[0], ldu, // U
                new double[0], ldvt, // Vt
                worksize, -1, new int[8 * min], info);

        int lwork;
        if (info.val == 0) {
            lwork = (int) worksize[0];
        } else {
            if (jobz == JobSVD.None) { // N
                lwork = 3 * min + Math.max(max, 7 * min);
            } else if (jobz == JobSVD.Part) { // S
                lwork = 3 * min * min + Math.max(max, 4 * min * min + 4 * min);
            } else if (jobz == JobSVD.Overwrite) { // O
                lwork = 3 * min * min + Math.max(max, 5 * min * min + 4 * min);
            } else { // All //A
                lwork = 3 * min * min + Math.max(max, 4 * min * min + 4 * min);
            }
        }

        lwork = Math.max(lwork, 1);

        return lwork;
    }

    /**
     * Convenience method for computing a full SVD
     * 
     * @param A
     *            Matrix to decompose, not modified
     * @return Newly allocated factorization
     * @throws NotConvergedException
     */
    public static SVD factorize(Matrix A) throws NotConvergedException {
        return new SVD(A.numRows(), A.numColumns()).factor(new DenseMatrix(A));
    }

    /**
     * Computes an SVD
     * 
     * @param A
     *            Matrix to decompose. Size must conform, and it will be
     *            overwritten on return. Pass a copy to avoid this
     * @return The current decomposition
     * @throws NotConvergedException
     */
    public SVD factor(DenseMatrix A) throws NotConvergedException {
        if (A.numRows() != m)
            throw new IllegalArgumentException("A.numRows() != m");
        else if (A.numColumns() != n)
            throw new IllegalArgumentException("A.numColumns() != n");

        intW info = new intW(0);

        if (job2 == null) {
            LAPACK.getInstance().dgesdd(job1.netlib(), m, n, A.getData(), lda,
                    S, U.getData(), ldu, Vt.getData(), ldvt, work, work.length,
                    iwork, info);

            if (job1 == JobSVD.Overwrite) {
                if (A.numRows >= A.numColumns) {
                    U = A;
                } else {
                    Vt = A;
                }
            }

        } else {
            LAPACK.getInstance().dgesvd(job1.netlib(), job2.netlib(), m, n,
                    A.getData(), lda, // A
                    S, // S
                    U.getData(), ldu, // U
                    Vt.getData(), ldvt, // Vt
                    work, work.length, info);

            if (job1 == JobSVD.Overwrite) {
                U = A;
            }

            if (job2 == JobSVD.Overwrite) {
                Vt = A;
            }
        }

        if (info.val > 0)
            throw new NotConvergedException(
                    NotConvergedException.Reason.Iterations);
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return this;
    }

    /**
     * True if singular vectors are stored
     */
    public boolean hasSingularVectors() {
        return job2 == null ? job1 != JobSVD.None : job1 != JobSVD.None
                || job2 != JobSVD.None;
    }

    /**
     * True if left singular vectors (U) are stored
     */
    public boolean hasLeftSingularVectors() {
        return job1 != JobSVD.None;
    }

    /**
     * True if right singular vectors (Vt) are stored
     */
    public boolean hasRightSingularVectors() {
        return job2 == null ? job1 != JobSVD.None : job2 != JobSVD.None;
    }

    /**
     * Returns the left singular vectors, column-wise. Not available for partial
     * decompositions
     * 
     * @return Matrix of size m*m
     */
    public DenseMatrix getU() {
        return U;
    }

    /**
     * Returns the right singular vectors, row-wise. Not available for partial
     * decompositions
     * 
     * @return Matrix of size n*n
     */
    public DenseMatrix getVt() {
        return Vt;
    }

    /**
     * Returns the singular values (stored in descending order)
     * 
     * @return Array of size min(m,n)
     */
    public double[] getS() {
        return S;
    }

}
