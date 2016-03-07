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

import org.junit.Test;

/**
 * QL test
 */
public class QLTest extends OrthogonalTestAbstract {

    @Test
    public void testStaticFactorize() {
        assertEqualsQL(A, QL.factorize(A));
    }

    @Test
    public void testRepeatStaticFactorize() {
        assertEqualsQL(A, QL.factorize(A));
        assertEqualsQL(A, QL.factorize(A));
    }

    @Test
    public void testFactor() {
        QL ql = new QL(A.numRows(), A.numColumns());
        assertEqualsQL(A, ql.factor(new DenseMatrix(A)));
    }

    @Test
    public void testRepeatFactor() {
        QL ql = new QL(A.numRows(), A.numColumns());
        ql.factor(new DenseMatrix(A));
        assertEqualsQL(A, ql);
        ql.factor(new DenseMatrix(A));
        assertEqualsQL(A, ql);
    }

    @Test
    public void testStaticFactorizeNonSquare() {
        assertEqualsQL(Ar, QL.factorize(Ar));
    }

    @Test
    public void testRepeatStaticFactorizeNonSquare() {
        assertEqualsQL(Ar, QL.factorize(Ar));
        assertEqualsQL(Ar, QL.factorize(Ar));
    }

    @Test
    public void testFactorNonSquare() {
        QL ql = new QL(Ar.numRows(), Ar.numColumns());
        assertEqualsQL(Ar, ql.factor(new DenseMatrix(Ar)));
    }

    @Test
    public void testRepeatFactorNonSquare() {
        QL ql = new QL(Ar.numRows(), Ar.numColumns());
        ql.factor(new DenseMatrix(Ar));
        assertEqualsQL(Ar, ql);
        ql.factor(new DenseMatrix(Ar));
        assertEqualsQL(Ar, ql);
    }

    private void assertEqualsQL(Matrix A, QL ql) {
        assertMatrixEquals(A, ql.getQ().mult(ql.getL(), A.copy().zero()));
    }

}
