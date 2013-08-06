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

/**
 * Base class for the orthogonal matrix decompositions (QR, RQ, LQ, and QL)
 */
abstract class OrthogonalComputer {

    /**
     * The orthogonal matrix
     */
    final DenseMatrix Q;

    /**
     * Lower triangular factor. May not be present
     */
    final LowerTriangDenseMatrix L;

    /**
     * Upper triangular factor. May not be present
     */
    final UpperTriangDenseMatrix R;

    /**
     * Factorisation sizes
     */
    final int m, n, k;

    /**
     * Work arrays
     */
    double[] work, workGen;

    /**
     * Scales for the reflectors
     */
    final double[] tau;

    /**
     * Constructor for OrthogonalComputer
     * 
     * @param m
     *            Number of rows
     * @param n
     *            Number of columns
     * @param upper
     *            True for storing an upper triangular factor, false for a lower
     *            triangular factor
     */
    OrthogonalComputer(int m, int n, boolean upper) {
        this.m = m;
        this.n = n;
        this.k = Math.min(m, n);

        tau = new double[k];

        Q = new DenseMatrix(m, n);
        if (upper) {
            R = new UpperTriangDenseMatrix(Math.min(m, n));
            L = null;
        } else {
            L = new LowerTriangDenseMatrix(Math.min(m, n));
            R = null;
        }
    }

    /**
     * Computes an orthogonal decomposition
     * 
     * @param A
     *            Matrix to decompose. Overwritten on exit. Pass a copy to avoid
     *            this
     * @return The current decomposition
     */
    public abstract OrthogonalComputer factor(DenseMatrix A);

    /**
     * Returns the orthogonal part of the factorization
     */
    public DenseMatrix getQ() {
        return Q;
    }

}
