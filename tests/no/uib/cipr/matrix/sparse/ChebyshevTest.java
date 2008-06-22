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

import no.uib.cipr.matrix.SymmDenseEVD;
import no.uib.cipr.matrix.sparse.Chebyshev;

/**
 * Test of Chebyshev solver
 */
public class ChebyshevTest extends SPDIterativeSolverTestAbstract {

    public ChebyshevTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void createSolver() throws Exception {
        // Get the extremal eigenvalues
        SymmDenseEVD evd = SymmDenseEVD.factorize(A);
        double[] eigs = evd.getEigenvalues();

        double eigmin = 1, eigmax = 1;
        if (eigs.length > 0) {
            eigmin = eigs[0];
            eigmax = eigs[eigs.length - 1];
        }

        solver = new Chebyshev(x, eigmin, eigmax);
    }

}
