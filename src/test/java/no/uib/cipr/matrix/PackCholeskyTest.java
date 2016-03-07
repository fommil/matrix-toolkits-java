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
 * Tests the packed Cholesky decomposition
 */
public class PackCholeskyTest {

    private LowerSPDPackMatrix L;

    private UpperSPDPackMatrix U;

    private DenseMatrix I;

    private final int max = 50;

    @Before
    public void setUp() throws Exception {
        int n = Utilities.getInt(1, max);

        L = new LowerSPDPackMatrix(n);
        Utilities.lowerPopulate(L);
        Utilities.addDiagonal(L, 1);
        while (!Utilities.spd(L))
            Utilities.addDiagonal(L, 1);

        U = new UpperSPDPackMatrix(n);
        Utilities.upperPopulate(U);
        Utilities.addDiagonal(U, 1);
        while (!Utilities.spd(U))
            Utilities.addDiagonal(U, 1);

        I = Matrices.identity(n);
    }

    @After
    public void tearDown() throws Exception {
        L = null;
        U = null;
        I = null;
    }

    @Test
    public void testLowerPackCholesky() {
        int n = L.numRows();

        PackCholesky c = new PackCholesky(n, false);
        c.factor(L.copy());

        c.solve(I);

        Matrix J = I.mult(L, new DenseMatrix(n, n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j)
                    assertEquals(J.get(i, j), 0, 1e-10);
                else
                    assertEquals(J.get(i, j), 1, 1e-10);
    }

    @Test
    public void testUpperPackCholesky() {
        int n = U.numRows();

        PackCholesky c = new PackCholesky(n, true);
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

    @Test
    public void testLowerPackCholeskyrcond() {
        int n = L.numRows();

        PackCholesky c = new PackCholesky(n, false);
        c.factor(L.copy());

        c.rcond(L);
    }

    @Test
    public void testUpperPackCholeskyrcond() {
        int n = U.numRows();

        PackCholesky c = new PackCholesky(n, true);
        c.factor(U.copy());

        c.rcond(U);
    }
}
