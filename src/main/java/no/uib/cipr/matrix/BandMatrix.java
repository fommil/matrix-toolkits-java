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

import java.util.Arrays;

import com.github.fommil.netlib.BLAS;
import com.github.fommil.netlib.LAPACK;
import org.netlib.util.intW;


/**
 * Banded matrix. The banded matrix is a useful sparse structure for many kinds
 * of direct computations, however it should only be used if the band is
 * sufficiently narrow as wide bands actually wastes both memory and compute
 * time. The matrix
 * <p>
 * <table border="1">
 * <tr>
 * <td>a<sub>11</sub></td>
 * <td>a<sub>12</sub></td>
 * <td>&nbsp;</td>
 * <td>&nbsp;</td>
 * <td>&nbsp;</td>
 * </tr>
 * <tr>
 * <td>a<sub>21</sub></td>
 * <td>a<sub>22</sub></td>
 * <td>a<sub>23</sub></td>
 * <td>&nbsp;</td>
 * <td>&nbsp;</td>
 * </tr>
 * <tr>
 * <td>a<sub>31</sub></td>
 * <td>a<sub>32</sub></td>
 * <td>a<sub>33</sub></td>
 * <td>a<sub>34</sub></td>
 * <td>&nbsp;</td>
 * </tr>
 * <tr>
 * <td>&nbsp;</td>
 * <td>a<sub>42</sub></td>
 * <td>a<sub>43</sub></td>
 * <td>a<sub>44</sub></td>
 * <td>a<sub>45</sub></td>
 * </tr>
 * <tr>
 * <td>&nbsp;</td>
 * <td>&nbsp;</td>
 * <td>a<sub>53</sub></td>
 * <td>a<sub>54</sub></td>
 * <td>a<sub>55</sub></td>
 * </tr>
 * </table>
 * </p>
 * <p>
 * has two lower diagonals and one upper diagonal. It will be stored in the
 * array
 * </p>
 * <p>
 * <table border="1">
 * <tr>
 * <td>&nbsp;</td>
 * <td>a<sub>11</sub></td>
 * <td>a<sub>21</sub></td>
 * <td>a<sub>31</sub></td>
 * <td>a<sub>21</sub></td>
 * <td>a<sub>22</sub></td>
 * <td>a<sub>32</sub></td>
 * <td>a<sub>42</sub></td>
 * <td>a<sub>23</sub></td>
 * <td>a<sub>33</sub></td>
 * <td>a<sub>43</sub></td>
 * <td>a<sub>53</sub></td>
 * <td>a<sub>34</sub></td>
 * <td>a<sub>44</sub></td>
 * <td>a<sub>54</sub></td>
 * <td>&nbsp;</td>
 * <td>a<sub>45</sub></td>
 * <td>a<sub>55</sub></td>
 * <td>&nbsp;</td>
 * <td>&nbsp;</td>
 * </tr>
 * </table>
 * </p>
 * <p>
 * Empty cells are allocated, but never referenced.
 * </p>
 */
public class BandMatrix extends AbstractBandMatrix {

    /**
     * Constructor for BandMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     * @param kl
     *            Number of bands above the main diagonal (superdiagonals)
     * @param ku
     *            Number of bands below the main diagonal (subdiagonals)
     */
    public BandMatrix(int n, int kl, int ku) {
        super(n, kl, ku);
    }

    /**
     * Constructor for BandMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the parts of <code>A</code>
     *            that lie within the allocated band are copied over, the rest
     *            is ignored
     * @param kl
     *            Number of bands above the main diagonal (superdiagonals)
     * @param ku
     *            Number of bands below the main diagonal (subdiagonals)
     */
    public BandMatrix(Matrix A, int kl, int ku) {
        super(A, kl, ku);
    }

    /**
     * Constructor for BandMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only the parts of <code>A</code>
     *            that lie within the allocated band are copied over, the rest
     *            is ignored
     * @param kl
     *            Number of bands above the main diagonal (superdiagonals)
     * @param ku
     *            Number of bands below the main diagonal (subdiagonals)
     * @param deep
     *            True for a deep copy. For shallow copies, <code>A</code>
     *            must be a banded matrix
     */
    public BandMatrix(Matrix A, int kl, int ku, boolean deep) {
        super(A, kl, ku, deep);
    }

    @Override
    public BandMatrix copy() {
        return new BandMatrix(this, kl, ku);
    }

    @Override
    public Matrix zero() {
        Arrays.fill(data, 0);
        return this;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.multAdd(alpha, x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dgbmv(Transpose.NoTranspose.netlib(), numRows, numColumns, kl,
                ku, alpha, data, kl + ku + 1, xd, 1, 1, yd, 1);

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.transMultAdd(alpha, x, y);

        checkTransMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dgbmv(Transpose.Transpose.netlib(), numRows, numColumns, kl, ku,
                alpha, data, kl + ku + 1, xd, 1, 1, yd, 1);

        return y;
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        if (!(X instanceof DenseMatrix))
            throw new UnsupportedOperationException("X must be a DenseMatrix");

        checkSolve(B, X);

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        // Allocate factorization matrix. The factorization matrix will be
        // large enough to accomodate any pivots
        BandMatrix Af = new BandMatrix(this, kl, ku + kl);
        int[] ipiv = new int[numRows];

        intW info = new intW(0);
        LAPACK.getInstance().dgbsv(numRows, kl, ku, X.numColumns(),
                Af.getData(), Matrices.ld(2 * kl + ku + 1), ipiv, Xd, Matrices.ld(numRows), info);

        if (info.val > 0)
            throw new MatrixSingularException();
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return X;
    }

    @Override
    public Vector solve(Vector b, Vector x) {
        DenseMatrix B = new DenseMatrix(b, false), X = new DenseMatrix(x, false);
        solve(B, X);
        return x;
    }

    @Override
    public Matrix transpose() {
        checkTranspose();

        if (kl != ku)
            throw new IllegalArgumentException("kl != ku");

        for (int j = 0; j < numColumns; ++j)
            for (int i = j + 1; i < Math.min(j + kl + 1, numRows); ++i) {
                double value = get(i, j);
                set(i, j, get(j, i));
                set(j, i, value);
            }

        return this;
    }

}
