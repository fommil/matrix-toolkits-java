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
 * Symmetric eigenvalue decomposition
 */
abstract class SymmEVD {

    /**
     * Size of the matrix
     */
    final int n;

    /**
     * The eigenvalues
     */
    final double[] w;

    /**
     * The eigenvectors stored columnwise
     */
    final DenseMatrix Z;

    /**
     * Job to do
     */
    final JobEig job;

    /**
     * Allocates storage for an eigenvalue computation
     * 
     * @param n
     *            Size of the matrix
     * @param vectors
     *            True to compute the eigenvectors, false for just the
     *            eigenvalues
     */
    public SymmEVD(int n, boolean vectors) {
        this.n = n;
        w = new double[n];
        job = vectors ? JobEig.All : JobEig.Eigenvalues;

        if (vectors)
            Z = new DenseMatrix(n, n);
        else
            Z = null;
    }

    /**
     * Allocates storage for an eigenvalue computation. Includes eigenvectors
     * 
     * @param n
     *            Size of the matrix
     */
    public SymmEVD(int n) {
        this(n, true);
    }

    /**
     * Gets the eigenvalues (stored in ascending order)
     */
    public double[] getEigenvalues() {
        return w;
    }

    /**
     * Gets the eigenvectors, if available
     */
    public DenseMatrix getEigenvectors() {
        return Z;
    }

    /**
     * True if the eigenvectors have been computed
     */
    public boolean hasEigenvectors() {
        return Z != null;
    }

}
