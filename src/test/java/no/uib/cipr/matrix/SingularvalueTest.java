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

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Test the singular value solver
 */
public class SingularvalueTest {

    /**
     * Matrix to decompose
     */
    private DenseMatrix A;

    /**
     * Maximum matrix size, to avoid too slow tests
     */
    private final int max = 100;

    @Before
    public void setUp() throws Exception {
        int n = Utilities.getInt(1, max);
        A = new DenseMatrix(n, n);
    }

    @After
    public void tearDown() throws Exception {
        A = null;
    }

    @Test
    public void testStaticFactorize() throws NotConvergedException {
        assertEqualsSVD(A, SVD.factorize(A));
    }

    @Test
    public void testFactor() throws NotConvergedException {
        SVD svd = new SVD(A.numRows(), A.numColumns());
        assertEqualsSVD(A, svd.factor(A.copy()));
    }

    private void assertEqualsSVD(Matrix A, SVD svd) {
        TridiagMatrix S = new TridiagMatrix(svd.getS().length);
        System.arraycopy(svd.getS(), 0, S.getDiagonal(), 0, svd.getS().length);
        DenseMatrix U = svd.getU();
        DenseMatrix Vt = svd.getVt();

        // Compute U*S*Vt
        Matrix s = U.mult(
                S.mult(Vt, new DenseMatrix(S.numRows(), Vt.numColumns())),
                new DenseMatrix(A.numRows(), A.numColumns()));

        // Check that A=U*S*Vt
        for (int i = 0; i < A.numRows(); ++i)
            for (int j = 0; j < A.numColumns(); ++j)
                assertEquals(A.get(i, j), s.get(i, j), 1e-12);
    }

}
