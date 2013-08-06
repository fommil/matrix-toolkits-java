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
 * Upper triangular banded matrix. The same storage as
 * {@link no.uib.cipr.matrix.BandMatrix BandMatrix}, but without subdiagonals.
 */
public class UpperTriangBandMatrix extends AbstractTriangBandMatrix {

    /**
     * Constructor for UpperTriangBandMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     * @param kd
     *            Number of bands above the main diagonal (superdiagonals)
     */
    public UpperTriangBandMatrix(int n, int kd) {
        super(n, 0, kd, UpLo.Upper, Diag.NonUnit);
    }

    /**
     * Constructor for UpperTriangBandMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the parts of <code>A</code>
     *            that lie within the allocated band are copied over, the rest
     *            is ignored
     * @param kd
     *            Number of bands above the main diagonal (superdiagonals)
     */
    public UpperTriangBandMatrix(Matrix A, int kd) {
        this(A, kd, true);
    }

    /**
     * Constructor for UpperTriangBandMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the parts of <code>A</code>
     *            that lie within the allocated band are copied over, the rest
     *            is ignored
     * @param kd
     *            Number of bands above the main diagonal (superdiagonals)
     * @param deep
     *            True for a deep copy. For shallow copies, <code>A</code>
     *            must be a banded matrix
     */
    public UpperTriangBandMatrix(Matrix A, int kd, boolean deep) {
        super(A, 0, kd, deep, UpLo.Upper, Diag.NonUnit);
    }

    /**
     * Constructor for UpperTriangBandMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     * @param kd
     *            Number of bands above the main diagonal (superdiagonals)
     */
    UpperTriangBandMatrix(int n, int kd, Diag diag) {
        super(n, 0, kd, UpLo.Upper, diag);
    }

    /**
     * Constructor for UpperTriangBandMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the parts of <code>A</code>
     *            that lie within the allocated band are copied over, the rest
     *            is ignored
     * @param kd
     *            Number of bands above the main diagonal (superdiagonals)
     * @param deep
     *            True for a deep copy. For shallow copies, <code>A</code>
     *            must be a banded matrix
     */
    UpperTriangBandMatrix(Matrix A, int kd, boolean deep, Diag diag) {
        super(A, 0, kd, deep, UpLo.Upper, diag);
    }

    @Override
    public UpperTriangBandMatrix copy() {
        return new UpperTriangBandMatrix(this, ku);
    }

}
