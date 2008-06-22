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

import no.uib.cipr.matrix.LowerSymmDenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SymmDenseEVD;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;

/**
 * Test of the symmetric, dense eigenvalue solver
 */
public class SymmDenseEigenvalueTest extends SymmEigenvalueTestAbstract {

    private LowerSymmDenseMatrix L;

    private UpperSymmDenseMatrix U;

    public SymmDenseEigenvalueTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        L = new LowerSymmDenseMatrix(A);
        U = new UpperSymmDenseMatrix(A);
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
        L = null;
        U = null;
    }

    public void testLowerStaticFactorize() throws NotConvergedException {
        SymmDenseEVD evd = SymmDenseEVD.factorize(L);
        assertEquals(L, evd.getEigenvalues(), evd.getEigenvectors());
    }

    public void testUpperStaticFactorize() throws NotConvergedException {
        SymmDenseEVD evd = SymmDenseEVD.factorize(U);
        assertEquals(U, evd.getEigenvalues(), evd.getEigenvectors());
    }

    public void testLowerFactor() throws NotConvergedException {
        SymmDenseEVD evd = new SymmDenseEVD(A.numRows(), false);
        evd.factor(L.copy());
        assertEquals(L, evd.getEigenvalues(), evd.getEigenvectors());
    }

    public void testUpperFactor() throws NotConvergedException {
        SymmDenseEVD evd = new SymmDenseEVD(A.numRows(), true);
        evd.factor(U.copy());
        assertEquals(U, evd.getEigenvalues(), evd.getEigenvectors());
    }

    public void testLowerRepeatFactor() throws NotConvergedException {
        SymmDenseEVD evd = new SymmDenseEVD(A.numRows(), false);
        evd.factor(L.copy());
        assertEquals(L, evd.getEigenvalues(), evd.getEigenvectors());
        evd.factor(L.copy());
        assertEquals(L, evd.getEigenvalues(), evd.getEigenvectors());
    }

    public void testUpperRepeatFactor() throws NotConvergedException {
        SymmDenseEVD evd = new SymmDenseEVD(A.numRows(), true);
        evd.factor(U.copy());
        assertEquals(U, evd.getEigenvalues(), evd.getEigenvectors());
        evd.factor(U.copy());
        assertEquals(U, evd.getEigenvalues(), evd.getEigenvectors());
    }

}
