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
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import no.uib.cipr.matrix.Utilities;

/**
 * Test of iterative solvers for SPD matrices
 */
public abstract class SPDIterativeSolverTestAbstract extends IterativeSolverTestAbstract {

    public SPDIterativeSolverTestAbstract(String arg0) {
        super(arg0);
    }

    @Override
    protected void createMatrix() throws Exception {
        // Create a symmetrical matrix
        int n = Utilities.getInt(1, max);
        int b = Utilities.getInt(Math.min(bmax, n));
        A = new FlexCompRowMatrix(n, n);
        Utilities.symmetryPopulate(A, b);

        // Need positive eigenvalues
        addDiagonal(A, shift);
        SymmDenseEVD evd = SymmDenseEVD.factorize(A);
        while (n > 0 && evd.getEigenvalues()[0] <= 0) {
            addDiagonal(A, shift);
            evd = SymmDenseEVD.factorize(A);
        }
    }

}
