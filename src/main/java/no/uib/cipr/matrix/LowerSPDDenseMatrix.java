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
 * Lower symmetrical positive definite dense matrix. Same layout as
 * {@link no.uib.cipr.matrix.LowerSymmDenseMatrix LowerSymmDenseMatrix}. This
 * class does not enforce the SPD property, but serves as a tag so that more
 * efficient algorithms can be used in the solvers.
 */
public class LowerSPDDenseMatrix extends LowerSymmDenseMatrix {

    /**
     * Constructor for LowerSPDDenseMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public LowerSPDDenseMatrix(int n) {
        super(n);
    }

    /**
     * Constructor for LowerSPDDenseMatrix
     * 
     * @param A
     *            Matrix to copy. It must be a square matrix, and only the lower
     *            triangular part is copied
     */
    public LowerSPDDenseMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for LowerSPDDenseMatrix
     * 
     * @param A
     *            Matrix to copy. It must be a square matrix, and only the lower
     *            triangular part is copied
     * @param deep
     *            False for a shallow copy, else it'll be a deep copy. For
     *            shallow copies, <code>A</code> must be a dense matrix
     */
    public LowerSPDDenseMatrix(Matrix A, boolean deep) {
        super(A, deep);
    }

    @Override
    public LowerSPDDenseMatrix copy() {
        return new LowerSPDDenseMatrix(this);
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        return SPDsolve(B, X);
    }

}
