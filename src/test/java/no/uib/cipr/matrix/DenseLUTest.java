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
 * Tests the dense LU decomposition
 */
public class DenseLUTest {

    /**
     * Matrix to decompose
     */
    private DenseMatrix A;

    private DenseMatrix I;

    private final int max = 100;

    @Before
    public void setUp() throws Exception {
        int n = Utilities.getInt(1, max);
        A = new DenseMatrix(n, n);
        Utilities.populate(A);
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
    public void testDenseLU() {
        int n = A.numRows();
        DenseLU lu = new DenseLU(n, n);
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
    public void testDenseLUtranspose() {
        int n = A.numRows();
        DenseLU lu = new DenseLU(n, n);
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
    public void testDenseLUrcond() {
        int n = A.numRows();
        DenseLU lu = new DenseLU(n, n);
        lu.factor(A.copy());

        lu.rcond(A, Matrix.Norm.One);
        lu.rcond(A, Matrix.Norm.Infinity);
    }

    @Test
    public void testDensePLU() {
        Matrix m = new DenseMatrix(new double[][]{{2, -1, -2}, {-4, 6, 3},
                {-4, -2, 8}});
        DenseLU dlu = DenseLU.factorize(m);

        Matrix p = dlu.getP();
        Matrix l = dlu.getL();
        Matrix u = dlu.getU();

        Matrix lu = l.mult(u, new DenseMatrix(3, 3));
        Matrix x = p.mult(lu, new DenseMatrix(3, 3));

        MatrixTestAbstract.assertMatrixEquals(m, x);
    }

}
