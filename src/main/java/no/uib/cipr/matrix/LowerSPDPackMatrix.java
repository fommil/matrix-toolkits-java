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
 * Lower symmetrical positive definite packed matrix. Same layout as
 * {@link no.uib.cipr.matrix.LowerSymmPackMatrix LowerSymmPackMatrix}. This
 * class does not enforce the SPD property, but serves as a tag so that more
 * efficient algorithms can be used in the solvers.
 */
public class LowerSPDPackMatrix extends LowerSymmPackMatrix {

    /**
     * Constructor for LowerSPDPackMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public LowerSPDPackMatrix(int n) {
        super(n);
    }

    /**
     * Constructor for LowerSPDPackMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the entries of the relevant
     *            part are copied
     */
    public LowerSPDPackMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for LowerSPDPackMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the entries of the relevant
     *            part are copied
     * @param deep
     *            True if the copy is deep, else false (giving a shallow copy).
     *            For shallow copies, <code>A</code> must be a packed matrix
     */
    public LowerSPDPackMatrix(Matrix A, boolean deep) {
        super(A, deep);
    }

    @Override
    public LowerSPDPackMatrix copy() {
        return new LowerSPDPackMatrix(this);
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        return SPDsolve(B, X);
    }

}
