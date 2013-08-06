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

/**
 * Test of a matrix whose structure cannot change, but its numerical values can
 */
public abstract class StructImmutableMatrixTestAbstract extends MatrixTestAbstract {

    public StructImmutableMatrixTestAbstract(String arg0) {
        super(arg0);
    }

    @Override
    public void testMatrixAdd() {
        // Not supported
    }

    @Override
    public void testMatrixSet() {
        // Not supported
    }

    @Override
    public void testOneMatrixAdd() {
        // Not supported
    }

    @Override
    public void testOneMatrixSet() {
        // Not supported
    }

    @Override
    public void testRandomMatrixAdd() {
        // Not supported
    }

    @Override
    public void testRandomMatrixSet() {
        // Not supported
    }

    @Override
    public void testZeroMatrixAdd() {
        // Not supported
    }

    @Override
    public void testZeroMatrixSet() {
        // Not supported
    }

    @Override
    public void testVectorRank1() {
        // Not supported
    }

    @Override
    public void testVectorRank1Dense() {
        // Not supported
    }

    @Override
    public void testVectorRank2() {
        // Not supported
    }

    @Override
    public void testVectorRank2Dense() {
        // Not supported
    }

    @Override
    public void testMatrixRank1() {
        // Not supported
    }

    @Override
    public void testMatrixRank1Dense() {
        // Not supported
    }

    @Override
    public void testMatrixTransRank1Dense() {
        // Not supported
    }

    @Override
    public void testMatrixTransRank1() {
        // Not supported
    }

    @Override
    public void testMatrixRank2() {
        // Not supported
    }

    @Override
    public void testMatrixRank2Dense() {
        // Not supported
    }

    @Override
    public void testMatrixTransRank2() {
        // Not supported
    }

    @Override
    public void testMatrixTransRank2Dense() {
        // Not supported
    }

    @Override
    public void testTranspose() {
        // Not supported
    }

    @Override
    public void testTransposeInplace() {
        // Not supported
    }

}
