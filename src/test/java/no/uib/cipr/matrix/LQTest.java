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
 * LQ test
 */
public class LQTest extends OrthogonalTestAbstract {

    @Test
    public void testStaticFactorize() {
        assertLQEquals(A, LQ.factorize(A));
    }

    @Test
    public void testRepeatStaticFactorize() {
        assertLQEquals(A, LQ.factorize(A));
        assertLQEquals(A, LQ.factorize(A));
    }

    @Test
    public void testFactor() {
        LQ lq = new LQ(A.numRows(), A.numColumns());
        assertLQEquals(A, lq.factor(new DenseMatrix(A)));
    }

    @Test
    public void testRepeatFactor() {
        LQ lq = new LQ(A.numRows(), A.numColumns());
        lq.factor(new DenseMatrix(A));
        assertLQEquals(A, lq);
        lq.factor(new DenseMatrix(A));
        assertLQEquals(A, lq);
    }

    @Test
    public void testStaticFactorizeNonSquare() {
        assertLQEquals(Ac, LQ.factorize(Ac));
    }

    @Test
    public void testRepeatStaticFactorizeNonSquare() {
        assertLQEquals(Ac, LQ.factorize(Ac));
        assertLQEquals(Ac, LQ.factorize(Ac));
    }

    @Test
    public void testFactorNonSquare() {
        LQ lq = new LQ(Ac.numRows(), Ac.numColumns());
        assertLQEquals(Ac, lq.factor(new DenseMatrix(Ac)));
    }

    @Test
    public void testRepeatFactorNonSquare() {
        LQ lq = new LQ(Ac.numRows(), Ac.numColumns());
        lq.factor(new DenseMatrix(Ac));
        assertLQEquals(Ac, lq);
        lq.factor(new DenseMatrix(Ac));
        assertLQEquals(Ac, lq);
    }

    private void assertLQEquals(Matrix A, LQ lq) {
        assertMatrixEquals(A, lq.getL().mult(lq.getQ(), A.copy().zero()));
    }

}
