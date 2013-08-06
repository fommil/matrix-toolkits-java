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

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import junit.framework.TestCase;

/**
 * Tests the symmetric eigenvalue computers
 */
public abstract class SymmEigenvalueTestAbstract extends TestCase {

    /**
     * Initial work-matrix
     */
    protected Matrix A;

    /**
     * Maximum matrix size, to avoid too slow tests
     */
    private final int max = 100;

    public SymmEigenvalueTestAbstract(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
        int n = Utilities.getInt(max);
        A = Matrices.random(n, n);
    }

    @Override
    protected void tearDown() throws Exception {
        A = null;
    }

    protected void assertEquals(Matrix A, double[] w, DenseMatrix Z) {
        // A*X
        Matrix left = A.mult(Z, new DenseMatrix(A.numRows(), A.numColumns()));

        // lambda*X
        Matrix right = new DenseMatrix(Z);
        for (int i = 0; i < w.length; ++i)
            for (int j = 0; j < w.length; ++j)
                right.set(i, j, w[j] * right.get(i, j));

        // Check that A*X=lambda*X
        for (int i = 0; i < A.numRows(); ++i)
            for (int j = 0; j < A.numColumns(); ++j)
                assertEquals(left.get(i, j), right.get(i, j), 1e-12);
    }

}
