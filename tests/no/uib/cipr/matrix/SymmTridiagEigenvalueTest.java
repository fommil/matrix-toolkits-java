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

import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SymmTridiagEVD;
import no.uib.cipr.matrix.SymmTridiagMatrix;

/**
 * Test of the symmetric, tridiagonal eigenvalue solver
 */
public class SymmTridiagEigenvalueTest extends SymmEigenvalueTestAbstract {

    private SymmTridiagMatrix T;

    private final int max = 100;

    public SymmTridiagEigenvalueTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
        int n = Utilities.getInt(1, max);
        A = Matrices.random(n, n);
        T = new SymmTridiagMatrix(A);
    }

    public void testStaticFactorize() throws NotConvergedException {
        SymmTridiagEVD evd = SymmTridiagEVD.factorize(A);
        assertEquals(T, evd.getEigenvalues(), evd.getEigenvectors());
    }

    public void testFactor() throws NotConvergedException {
        SymmTridiagEVD evd = new SymmTridiagEVD(A.numRows());
        evd.factor(T.copy());
        assertEquals(T, evd.getEigenvalues(), evd.getEigenvectors());
    }

    public void testRepeatFactor() throws NotConvergedException {
        SymmTridiagEVD evd = new SymmTridiagEVD(A.numRows());
        evd.factor(T.copy());
        assertEquals(T, evd.getEigenvalues(), evd.getEigenvectors());
        evd.factor(T.copy());
        assertEquals(T, evd.getEigenvalues(), evd.getEigenvectors());
    }

}
