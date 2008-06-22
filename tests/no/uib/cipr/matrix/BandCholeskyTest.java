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

import junit.framework.TestCase;
import no.uib.cipr.matrix.BandCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.LowerSPDBandMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSPDBandMatrix;

/**
 * Tests the banded Cholesky decomposition
 */
public class BandCholeskyTest extends TestCase {

    private LowerSPDBandMatrix L;

    private UpperSPDBandMatrix U;

    private DenseMatrix I;

    private int kl, ku;

    private final int max = 100, bmax = 10;

    public BandCholeskyTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
        int n = Utilities.getInt(1, max);

        kl = Math.min(n, Utilities.getInt(bmax));
        ku = Math.min(n, Utilities.getInt(bmax));

        L = new LowerSPDBandMatrix(n, kl);
        Utilities.bandPopulate(L, kl, 0);
        Utilities.addDiagonal(L, 1);
        while (!Utilities.spd(L))
            Utilities.addDiagonal(L, 1);

        U = new UpperSPDBandMatrix(n, ku);
        Utilities.bandPopulate(U, 0, ku);
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

    public void testLowerBandCholesky() {
        int n = L.numRows();

        BandCholesky c = new BandCholesky(n, kl, false);
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

    public void testUpperBandCholesky() {
        int n = U.numRows();

        BandCholesky c = new BandCholesky(n, ku, true);
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

    public void testLowerBandCholeskyrcond() {
        int n = L.numRows();

        BandCholesky c = new BandCholesky(n, kl, false);
        c.factor(L.copy());

        c.rcond(L);
    }

    public void testUpperBandCholeskyrcond() {
        int n = U.numRows();

        BandCholesky c = new BandCholesky(n, ku, true);
        c.factor(U.copy());

        c.rcond(U);
    }
}
