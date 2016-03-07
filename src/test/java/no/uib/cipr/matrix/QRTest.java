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
 * QR test
 */
public class QRTest extends OrthogonalTestAbstract {

    @Test
    public void testStaticFactorize() {
        assertEqualsQR(A, QR.factorize(A));
    }

    @Test
    public void testRepeatStaticFactorize() {
        assertEqualsQR(A, QR.factorize(A));
        assertEqualsQR(A, QR.factorize(A));
    }

    @Test
    public void testFactor() {
        QR qr = new QR(A.numRows(), A.numColumns());
        assertEqualsQR(A, qr.factor(new DenseMatrix(A)));
    }

    @Test
    public void testRepeatFactor() {
        QR qr = new QR(A.numRows(), A.numColumns());
        qr.factor(new DenseMatrix(A));
        assertEqualsQR(A, qr);
        qr.factor(new DenseMatrix(A));
        assertEqualsQR(A, qr);
    }

    @Test
    public void testStaticFactorizeNonSquare() {
        assertEqualsQR(Ar, QR.factorize(Ar));
    }

    @Test
    public void testRepeatStaticFactorizeNonSquare() {
        assertEqualsQR(Ar, QR.factorize(Ar));
        assertEqualsQR(Ar, QR.factorize(Ar));
    }

    @Test
    public void testFactorNonSquare() {
        QR qr = new QR(Ar.numRows(), Ar.numColumns());
        assertEqualsQR(Ar, qr.factor(new DenseMatrix(Ar)));
    }

    @Test
    public void testRepeatFactorNonSquare() {
        QR qr = new QR(Ar.numRows(), Ar.numColumns());
        qr.factor(new DenseMatrix(Ar));
        assertEqualsQR(Ar, qr);
        qr.factor(new DenseMatrix(Ar));
        assertEqualsQR(Ar, qr);
    }

    private void assertEqualsQR(Matrix A, QR qr) {
        assertMatrixEquals(A, qr.getQ().mult(qr.getR(), A.copy().zero()));
    }

}
