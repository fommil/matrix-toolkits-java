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

import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import junit.framework.TestCase;

/**
 * Orthogonal matrix decomposition tests
 */
public abstract class OrthogonalTestAbstract extends TestCase {

    /**
     * Initial work-matrix, and non-square matrices
     */
    protected Matrix A, Ar, Ac;

    /**
     * Maximum matrix size, to avoid too slow tests
     */
    private final int max = 100;

    public OrthogonalTestAbstract(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
        int n = Utilities.getInt(1, max);
        int m = Utilities.getInt(1, n);
        A = Matrices.random(n, n);
        Ar = Matrices.random(n, m);
        Ac = Matrices.random(m, n);
    }

    @Override
    protected void tearDown() throws Exception {
        A = Ar = Ac = null;
    }

    protected void assertEquals(Matrix A, Matrix B) {
        for (int i = 0; i < A.numRows(); ++i)
            for (int j = 0; j < A.numColumns(); ++j)
                assertEquals(A.get(i, j), B.get(i, j), 1e-12);
    }

}
