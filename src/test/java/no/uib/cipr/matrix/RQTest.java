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
import no.uib.cipr.matrix.RQ;

/**
 * RQ test
 */
public class RQTest extends OrthogonalTestAbstract {

    public RQTest(String arg0) {
        super(arg0);
    }

    public void testStaticFactorize() {
        assertEquals(A, RQ.factorize(A));
    }

    public void testRepeatStaticFactorize() {
        assertEquals(A, RQ.factorize(A));
        assertEquals(A, RQ.factorize(A));
    }

    public void testFactor() {
        RQ c = new RQ(A.numRows(), A.numColumns());
        assertEquals(A, c.factor(new DenseMatrix(A)));
    }

    public void testRepeatFactor() {
        RQ rq = new RQ(A.numRows(), A.numColumns());
        rq.factor(new DenseMatrix(A));
        assertEquals(A, rq);
        rq.factor(new DenseMatrix(A));
        assertEquals(A, rq);
    }

    public void testStaticFactorizeNonSquare() {
        assertEquals(Ac, RQ.factorize(Ac));
    }

    public void testRepeatStaticFactorizeNonSquare() {
        assertEquals(Ac, RQ.factorize(Ac));
        assertEquals(Ac, RQ.factorize(Ac));
    }

    public void testFactorNonSquare() {
        RQ rq = new RQ(Ac.numRows(), Ac.numColumns());
        assertEquals(Ac, rq.factor(new DenseMatrix(Ac)));
    }

    public void testRepeatFactorNonSquare() {
        RQ rq = new RQ(Ac.numRows(), Ac.numColumns());
        rq.factor(new DenseMatrix(Ac));
        assertEquals(Ac, rq);
        rq.factor(new DenseMatrix(Ac));
        assertEquals(Ac, rq);
    }

    private void assertEquals(Matrix A, RQ rq) {
        assertEquals(A, rq.getR().mult(rq.getQ(), A.copy().zero()));
    }

}
