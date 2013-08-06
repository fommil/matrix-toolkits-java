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
import no.uib.cipr.matrix.QR;

/**
 * QR test
 */
public class QRTest extends OrthogonalTestAbstract {

    public QRTest(String arg0) {
        super(arg0);
    }

    public void testStaticFactorize() {
        assertEquals(A, QR.factorize(A));
    }

    public void testRepeatStaticFactorize() {
        assertEquals(A, QR.factorize(A));
        assertEquals(A, QR.factorize(A));
    }

    public void testFactor() {
        QR qr = new QR(A.numRows(), A.numColumns());
        assertEquals(A, qr.factor(new DenseMatrix(A)));
    }

    public void testRepeatFactor() {
        QR qr = new QR(A.numRows(), A.numColumns());
        qr.factor(new DenseMatrix(A));
        assertEquals(A, qr);
        qr.factor(new DenseMatrix(A));
        assertEquals(A, qr);
    }

    public void testStaticFactorizeNonSquare() {
        assertEquals(Ar, QR.factorize(Ar));
    }

    public void testRepeatStaticFactorizeNonSquare() {
        assertEquals(Ar, QR.factorize(Ar));
        assertEquals(Ar, QR.factorize(Ar));
    }

    public void testFactorNonSquare() {
        QR qr = new QR(Ar.numRows(), Ar.numColumns());
        assertEquals(Ar, qr.factor(new DenseMatrix(Ar)));
    }

    public void testRepeatFactorNonSquare() {
        QR qr = new QR(Ar.numRows(), Ar.numColumns());
        qr.factor(new DenseMatrix(Ar));
        assertEquals(Ar, qr);
        qr.factor(new DenseMatrix(Ar));
        assertEquals(Ar, qr);
    }

    private void assertEquals(Matrix A, QR qr) {
        assertEquals(A, qr.getQ().mult(qr.getR(), A.copy().zero()));
    }

}
