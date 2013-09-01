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

import com.github.fommil.netlib.LAPACK;
import org.netlib.util.intW;

/**
 * Symmetrical positive definite tridiagonal matrix. Same as
 * {@link no.uib.cipr.matrix.SymmTridiagMatrix SymmTridiagMatrix}, and is used
 * as a marker class to allow for more efficient solvers.
 */
public class SPDTridiagMatrix extends SymmTridiagMatrix {

    /**
     * Constructor for SPDTridiagMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public SPDTridiagMatrix(int n) {
        super(n);
    }

    /**
     * Constructor for SPDTridiagMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only main and the superdiagonal
     *            is copied over
     */
    public SPDTridiagMatrix(Matrix A) {
        super(A);
    }

    /**
     * Constructor for SPDTridiagMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only main and the superdiagonal
     *            is copied over
     * @param deep
     *            True for a deep copy. For shallow copies <code>A</code> must
     *            be a <code>SymmTridiagMatrix</code>
     */
    public SPDTridiagMatrix(Matrix A, boolean deep) {
        super(A, deep);
    }

    @Override
    public SPDTridiagMatrix copy() {
        return new SPDTridiagMatrix(this);
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        if (!(X instanceof DenseMatrix))
            throw new UnsupportedOperationException("X must be a DenseMatrix");

        checkSolve(B, X);

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        intW info = new intW(0);
        LAPACK.getInstance().dptsv(numRows, X.numColumns(),
                diag.clone(), offDiag.clone(), Xd, Matrices.ld(numRows), info);

        if (info.val > 0)
            throw new MatrixNotSPDException();
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return X;
    }

}
