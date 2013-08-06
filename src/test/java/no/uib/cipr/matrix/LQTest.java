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
import no.uib.cipr.matrix.LQ;
import no.uib.cipr.matrix.Matrix;

/**
 * LQ test
 */
public class LQTest extends OrthogonalTestAbstract {

    public LQTest(String arg0) {
        super(arg0);
    }

    public void testStaticFactorize() {
        assertEquals(A, LQ.factorize(A));
    }

    public void testRepeatStaticFactorize() {
        assertEquals(A, LQ.factorize(A));
        assertEquals(A, LQ.factorize(A));
    }

    public void testFactor() {
        LQ lq = new LQ(A.numRows(), A.numColumns());
        assertEquals(A, lq.factor(new DenseMatrix(A)));
    }

    public void testRepeatFactor() {
        LQ lq = new LQ(A.numRows(), A.numColumns());
        lq.factor(new DenseMatrix(A));
        assertEquals(A, lq);
        lq.factor(new DenseMatrix(A));
        assertEquals(A, lq);
    }

    public void testStaticFactorizeNonSquare() {
        assertEquals(Ac, LQ.factorize(Ac));
    }

    public void testRepeatStaticFactorizeNonSquare() {
        assertEquals(Ac, LQ.factorize(Ac));
        assertEquals(Ac, LQ.factorize(Ac));
    }

    public void testFactorNonSquare() {
        LQ lq = new LQ(Ac.numRows(), Ac.numColumns());
        assertEquals(Ac, lq.factor(new DenseMatrix(Ac)));
    }

    public void testRepeatFactorNonSquare() {
        LQ lq = new LQ(Ac.numRows(), Ac.numColumns());
        lq.factor(new DenseMatrix(Ac));
        assertEquals(Ac, lq);
        lq.factor(new DenseMatrix(Ac));
        assertEquals(Ac, lq);
    }

    private void assertEquals(Matrix A, LQ lq) {
        assertEquals(A, lq.getL().mult(lq.getQ(), A.copy().zero()));
    }

}
