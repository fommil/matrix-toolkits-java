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

import no.uib.cipr.matrix.StructImmutableMatrixTestAbstract;

/**
 * Test of sparse matrices whose sparsity structure is immutable
 */
public abstract class SparseStructImmutableMatrixTestAbstract extends
        StructImmutableMatrixTestAbstract {

    protected int bmax = 100, tmax = 10;

    public SparseStructImmutableMatrixTestAbstract(String arg0) {
        super(arg0);
    }

    @Override
    public void testMatrixSolve() {
        // Not applicable
    }

    @Override
    public void testTransMatrixSolve() {
        // Not applicable
    }

    @Override
    public void testTransVectorSolve() {
        // Not applicable
    }

    @Override
    public void testVectorSolve() {
        // Not applicable
    }

}
