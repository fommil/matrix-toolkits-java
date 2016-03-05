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

import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

/**
 * Test of the symmetric, tridiagonal eigenvalue solver
 */
public class SymmBandEigenvalueTest extends SymmEigenvalueTestAbstract {

    private LowerSymmBandMatrix L;

    private UpperSymmBandMatrix U;

    @Before
    @Override
    public void setUp() throws Exception {
        super.setUp();
        int kd = Utilities.getInt(1, A.numRows());
        L = new LowerSymmBandMatrix(A, kd);
        U = new UpperSymmBandMatrix(A, kd);
    }

    @After
    @Override
    public void tearDown() throws Exception {
        super.tearDown();
        L = null;
        U = null;
    }

    @Test
    public void testLowerStaticFactorize() throws NotConvergedException {
        SymmBandEVD evd = SymmBandEVD.factorize(L, L.numSubDiagonals());
        assertEqualsWZ(L, evd.getEigenvalues(), evd.getEigenvectors());
    }

    @Test
    public void testUpperStaticFactorize() throws NotConvergedException {
        SymmBandEVD evd = SymmBandEVD.factorize(U, U.numSuperDiagonals());
        assertEqualsWZ(U, evd.getEigenvalues(), evd.getEigenvectors());
    }

    @Test
    public void testLowerFactor() throws NotConvergedException {
        SymmBandEVD evd = new SymmBandEVD(A.numRows(), false);
        evd.factor(L.copy());
        assertEqualsWZ(L, evd.getEigenvalues(), evd.getEigenvectors());
    }

    @Test
    public void testUpperFactor() throws NotConvergedException {
        SymmBandEVD evd = new SymmBandEVD(A.numRows(), true);
        evd.factor(U.copy());
        assertEqualsWZ(U, evd.getEigenvalues(), evd.getEigenvectors());
    }

    @Ignore
    @Test
    public void testLowerRepeatFactor() throws NotConvergedException {
        SymmBandEVD evd = new SymmBandEVD(A.numRows(), false);
        evd.factor(L.copy());
        assertEqualsWZ(L, evd.getEigenvalues(), evd.getEigenvectors());
        evd.factor(L.copy());
        assertEqualsWZ(L, evd.getEigenvalues(), evd.getEigenvectors());
    }

    @Ignore
    @Test
    public void testUpperRepeatFactor() throws NotConvergedException {
        SymmBandEVD evd = new SymmBandEVD(A.numRows(), true);
        evd.factor(U.copy());
        assertEqualsWZ(U, evd.getEigenvalues(), evd.getEigenvectors());
        evd.factor(U.copy());
        assertEqualsWZ(U, evd.getEigenvalues(), evd.getEigenvectors());
    }

}
