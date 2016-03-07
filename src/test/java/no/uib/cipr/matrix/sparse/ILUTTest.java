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

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Test of ILU with thresholding
 */
public class ILUTTest extends IncompleteFactorizationTestAbstract {

    @Override
    void testFactorization(Matrix A, Vector x) {
        Vector b = A.mult(x, x.copy());

        ILUT ilut = new ILUT(new FlexCompRowMatrix(A));
        ilut.setMatrix(A);
        ilut.apply(b, x);

        Vector r = A.multAdd(-1, x, b.copy());

        assertEquals(0, r.norm(Vector.Norm.TwoRobust), 1e-5);
    }

    /**
     * Test for diagInd values exceeding row data array length after dropTol is
     * applied causing ArrayOutOfBounds Error
     */
    @Test
    public void testILUTDropOutOfBoundsError() {
        int matrixSize = 100;
        FlexCompRowMatrix triDiagMatix = new FlexCompRowMatrix(matrixSize,
                matrixSize);

        double dropTol = 1e-6;
        int p = 50;

        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {
                if (i == j) {
                    triDiagMatix.set(i, j, 1);
                } else if (Math.abs(i - j) == 1) {
                    triDiagMatix.set(i, j, dropTol * 1e-2);
                }
            }
        }

        ILUT ilut = new ILUT(new FlexCompRowMatrix(triDiagMatix), dropTol, p);
        ilut.setMatrix(triDiagMatix);
    }

}
