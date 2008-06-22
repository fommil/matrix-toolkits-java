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

import no.uib.cipr.matrix.sparse.CompDiagMatrix;
import no.uib.cipr.matrix.Utilities;

/**
 * Test of CompDiagMatrix
 */
public class CompDiagMatrixTest extends SparseMatrixTestAbstract {

    public CompDiagMatrixTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void createPrimary() throws Exception {
        int n = Utilities.getInt(1, max);
        int m = Utilities.getInt(1, max);
        int b = Utilities.getInt(Math.min(bmax, m));
        A = new CompDiagMatrix(n, m);
        Ad = Utilities.rowPopulate(A, b);
    }

}