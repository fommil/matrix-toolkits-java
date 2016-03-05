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
 * Test of symmetrical matrices
 */
public abstract class SymmetricMatrixTestAbstract extends MatrixTestAbstract {

    @Test
    @Override
    public void testMatrixAdd() {
        // Not supported
    }

    @Test
    @Override
    public void testMatrixSet() {
        // Not supported
    }

    @Test
    @Override
    public void testOneMatrixAdd() {
        // Not supported
    }

    @Test
    @Override
    public void testOneMatrixSet() {
        // Not supported
    }

    @Test
    @Override
    public void testRandomMatrixAdd() {
        // Not supported
    }

    @Test
    @Override
    public void testRandomMatrixSet() {
        // Not supported
    }

    @Test
    @Override
    public void testVectorRank1() {
        if (A.isSquare()) {
            double alpha = Math.random();
            assertMatrixEquals(rank1(alpha, xdR, xdR), A.rank1(alpha, xR, xR));
        }
    }

    @Test
    @Override
    public void testVectorRank1Dense() {
        if (A.isSquare()) {
            double alpha = Math.random();
            assertMatrixEquals(rank1(alpha, xdR, xdR),
                    A.rank1(alpha, xDenseR, xDenseR));
        }
    }

}
