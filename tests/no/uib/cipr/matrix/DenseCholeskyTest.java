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

import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.LowerSPDDenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSPDDenseMatrix;
import junit.framework.TestCase;

/**
 * Tests the dense Cholesky decomposition
 */
public class DenseCholeskyTest extends TestCase {

    private LowerSPDDenseMatrix L;

    private UpperSPDDenseMatrix U;

    private DenseMatrix I;

    private final int max = 50;

    public DenseCholeskyTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
        int n = Utilities.getInt(1, max);

        L = new LowerSPDDenseMatrix(n);
        Utilities.lowerPopulate(L);
        Utilities.addDiagonal(L, 1);
        while (!Utilities.spd(L))
            Utilities.addDiagonal(L, 1);

        U = new UpperSPDDenseMatrix(n);
        Utilities.upperPopulate(U);
        Utilities.addDiagonal(U, 1);
        while (!Utilities.spd(U))
            Utilities.addDiagonal(U, 1);

        I = Matrices.identity(n);
    }

    @Override
    protected void tearDown() throws Exception {
        L = null;
        U = null;
        I = null;
    }

    public void testLowerDenseCholesky() {
        int n = L.numRows();

        DenseCholesky c = new DenseCholesky(n, false);
        c.factor(L.copy());

		assert I != null;
        c.solve(I);

        Matrix J = I.mult(L, new DenseMatrix(n, n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j)
                    assertEquals(J.get(i, j), 0, 1e-10);
                else
                    assertEquals(J.get(i, j), 1, 1e-10);
    }

    public void testUpperDenseCholesky() {
        int n = U.numRows();

        DenseCholesky c = new DenseCholesky(n, true);
        c.factor(U.copy());

        c.solve(I);

        Matrix J = I.mult(U, new DenseMatrix(n, n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j)
                    assertEquals(J.get(i, j), 0, 1e-10);
                else
                    assertEquals(J.get(i, j), 1, 1e-10);
    }

    public void testLowerDenseCholeskyrcond() {
        int n = L.numRows();

        DenseCholesky c = new DenseCholesky(n, false);
        c.factor(L.copy());

        c.rcond(L);
    }

    public void testUpperDenseCholeskyrcond() {
        int n = U.numRows();

        DenseCholesky c = new DenseCholesky(n, true);
        c.factor(U.copy());

        c.rcond(U);
    }
}
