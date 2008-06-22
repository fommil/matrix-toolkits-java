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

package no.uib.cipr.matrix.sparse;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Utilities;
import junit.framework.TestCase;

/**
 * Test of incomplete factorizations
 */
public abstract class IncompleteFactorizationTestAbstract extends TestCase {

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    public void testTriDiagonal() {
        int n = Utilities.getInt(1, 10);
        Matrix A = new DenseMatrix(n, n);
        Vector x = new DenseVector(n);

        for (int i = 0; i < n; ++i) {
            A.set(i, i, 2);
            x.set(i, 1);
        }
        for (int i = 0; i < n - 1; ++i) {
            A.set(i, i + 1, -1);
            A.set(i + 1, i, -1);
        }

        testFactorization(A, x);
    }

    public void testPentaDiagonal() {
        int n = Utilities.getInt(1, 10);
        Matrix A = new DenseMatrix(n, n);
        Vector x = new DenseVector(n);

        for (int i = 0; i < n; ++i) {
            A.set(i, i, 4);
            x.set(i, 1);
        }
        for (int i = 0; i < n - 1; ++i) {
            A.set(i, i + 1, -1);
            A.set(i + 1, i, -1);
        }
        for (int i = 0; i < n - 2; ++i) {
            A.set(i, i + 2, -1);
            A.set(i + 2, i, -1);
        }

        testFactorization(A, x);
    }

    abstract void testFactorization(Matrix A, Vector x);
}
