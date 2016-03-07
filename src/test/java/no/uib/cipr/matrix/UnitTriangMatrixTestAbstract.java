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
 * Test of unit, triangular matrices
 */
public abstract class UnitTriangMatrixTestAbstract
        extends
            TriangMatrixTestAbstract {

    @Test
    public void testAddDiagonal() {
        // Not applicable to unit triangular matrices
    }

    @Test
    public void testAddOneDiagonal() {
        // Not applicable to unit triangular matrices
    }

    @Test
    public void testAddZeroDiagonal() {
        // Not applicable to unit triangular matrices
    }

    @Test
    @Override
    public void testIteratorSet() {
        double alpha = Math.random();
        for (MatrixEntry e : A)
            if (e.row() != e.column())
                e.set(e.get() * alpha);
        assertMatrixEquals(Utilities.unitSet(scale(alpha)), A);
    }

    @Test
    @Override
    public void testIteratorSetGet() {
        // Not applicable to unit triangular matrices
    }

    @Test
    @Override
    public void testScale() {
        // Not applicable to unit triangular matrices
    }

    @Test
    @Override
    public void testZero() {
        // Not applicable to unit triangular matrices
    }

    @Test
    @Override
    public void testZeroScale() {
        // Not applicable to unit triangular matrices
    }

    /**
     * We can't zero, so we do without
     */
    @Test
    @Override
    public void testCopy() {
        Matrix Ac = A.copy();
        assertMatrixEquals(Ad, Ac);
    }

    @Test
    @Override
    public void testAdd() {
        double alpha = Math.random();
        for (MatrixEntry e : A)
            if (e.row() != e.column()) {
                A.add(e.row(), e.column(), alpha);
                A.add(e.row(), e.column(), -alpha);
            }
        assertMatrixEquals(Ad, A);
    }

}
