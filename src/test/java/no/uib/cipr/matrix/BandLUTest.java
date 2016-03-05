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
 * Tests the banded LU decomposition
 */
public class BandLUTest {

    private BandMatrix A;

    private DenseMatrix I;

    private int kl, ku;

    private final int max = 100, bmax = 10;

    @Before
    public void setUp() throws Exception {
        int n = Utilities.getInt(1, max);
        kl = Math.min(n, Utilities.getInt(bmax));
        ku = Math.min(n, Utilities.getInt(bmax));
        A = new BandMatrix(n, kl, kl + ku);
        Utilities.bandPopulate(A, kl, ku);
        Utilities.addDiagonal(A, 1);
        while (Utilities.singular(A))
            Utilities.addDiagonal(A, 1);

        I = Matrices.identity(n);
    }

    @After
    public void tearDown() throws Exception {
        A = null;
        I = null;
    }

    @Test
    public void testBandLU() {
        int n = A.numRows();

        BandLU lu = new BandLU(n, kl, ku);
        lu.factor(A.copy());

        lu.solve(I);

        Matrix J = I.mult(A, new DenseMatrix(n, n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j)
                    assertEquals(J.get(i, j), 0, 1e-10);
                else
                    assertEquals(J.get(i, j), 1, 1e-10);
    }

    @Test
    public void testBandLUtranspose() {
        int n = A.numRows();

        BandLU lu = new BandLU(n, kl, ku);
        lu.factor(A.copy());

        lu.transSolve(I);

        Matrix J = I.transAmult(A, new DenseMatrix(n, n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j)
                    assertEquals(J.get(i, j), 0, 1e-10);
                else
                    assertEquals(J.get(i, j), 1, 1e-10);
    }

    @Test
    public void testBandLUrcond() {
        int n = A.numRows();

        BandLU lu = new BandLU(n, kl, ku);
        lu.factor(A.copy());

        lu.rcond(A, Matrix.Norm.One);
        lu.rcond(A, Matrix.Norm.Infinity);
    }
}
