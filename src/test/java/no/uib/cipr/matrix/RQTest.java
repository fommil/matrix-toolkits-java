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
 * RQ test
 */
public class RQTest extends OrthogonalTestAbstract {

    @Test
    public void testStaticFactorize() {
        assertEqualsRQ(A, RQ.factorize(A));
    }

    @Test
    public void testRepeatStaticFactorize() {
        assertEqualsRQ(A, RQ.factorize(A));
        assertEqualsRQ(A, RQ.factorize(A));
    }

    @Test
    public void testFactor() {
        RQ c = new RQ(A.numRows(), A.numColumns());
        assertEqualsRQ(A, c.factor(new DenseMatrix(A)));
    }

    @Test
    public void testRepeatFactor() {
        RQ rq = new RQ(A.numRows(), A.numColumns());
        rq.factor(new DenseMatrix(A));
        assertEqualsRQ(A, rq);
        rq.factor(new DenseMatrix(A));
        assertEqualsRQ(A, rq);
    }

    @Test
    public void testStaticFactorizeNonSquare() {
        assertEqualsRQ(Ac, RQ.factorize(Ac));
    }

    @Test
    public void testRepeatStaticFactorizeNonSquare() {
        assertEqualsRQ(Ac, RQ.factorize(Ac));
        assertEqualsRQ(Ac, RQ.factorize(Ac));
    }

    @Test
    public void testFactorNonSquare() {
        RQ rq = new RQ(Ac.numRows(), Ac.numColumns());
        assertEqualsRQ(Ac, rq.factor(new DenseMatrix(Ac)));
    }

    @Test
    public void testRepeatFactorNonSquare() {
        RQ rq = new RQ(Ac.numRows(), Ac.numColumns());
        rq.factor(new DenseMatrix(Ac));
        assertEqualsRQ(Ac, rq);
        rq.factor(new DenseMatrix(Ac));
        assertEqualsRQ(Ac, rq);
    }

    private void assertEqualsRQ(Matrix A, RQ rq) {
        assertMatrixEquals(A, rq.getR().mult(rq.getQ(), A.copy().zero()));
    }

}
