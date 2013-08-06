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
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.QL;

/**
 * QL test
 */
public class QLTest extends OrthogonalTestAbstract {

    public QLTest(String arg0) {
        super(arg0);
    }

    public void testStaticFactorize() {
        assertEquals(A, QL.factorize(A));
    }

    public void testRepeatStaticFactorize() {
        assertEquals(A, QL.factorize(A));
        assertEquals(A, QL.factorize(A));
    }

    public void testFactor() {
        QL ql = new QL(A.numRows(), A.numColumns());
        assertEquals(A, ql.factor(new DenseMatrix(A)));
    }

    public void testRepeatFactor() {
        QL ql = new QL(A.numRows(), A.numColumns());
        ql.factor(new DenseMatrix(A));
        assertEquals(A, ql);
        ql.factor(new DenseMatrix(A));
        assertEquals(A, ql);
    }

    public void testStaticFactorizeNonSquare() {
        assertEquals(Ar, QL.factorize(Ar));
    }

    public void testRepeatStaticFactorizeNonSquare() {
        assertEquals(Ar, QL.factorize(Ar));
        assertEquals(Ar, QL.factorize(Ar));
    }

    public void testFactorNonSquare() {
        QL ql = new QL(Ar.numRows(), Ar.numColumns());
        assertEquals(Ar, ql.factor(new DenseMatrix(Ar)));
    }

    public void testRepeatFactorNonSquare() {
        QL ql = new QL(Ar.numRows(), Ar.numColumns());
        ql.factor(new DenseMatrix(Ar));
        assertEquals(Ar, ql);
        ql.factor(new DenseMatrix(Ar));
        assertEquals(Ar, ql);
    }

    private void assertEquals(Matrix A, QL ql) {
        assertEquals(A, ql.getQ().mult(ql.getL(), A.copy().zero()));
    }

}
