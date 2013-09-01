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

import java.io.IOException;

import no.uib.cipr.matrix.io.MatrixInfo;
import no.uib.cipr.matrix.io.MatrixSize;
import no.uib.cipr.matrix.io.MatrixVectorReader;

import com.github.fommil.netlib.BLAS;
import com.github.fommil.netlib.LAPACK;
import org.netlib.util.intW;

/**
 * Dense matrix. It is a good all-round matrix structure, with fast access and
 * efficient algebraic operations. The matrix
 * <p>
 * <table border="1">
 * <tr>
 * <td>a<sub>11</sub></td>
 * <td>a<sub>12</sub></td>
 * <td>a<sub>13</sub></td>
 * <td>a<sub>14</sub></td>
 * </tr>
 * <tr>
 * <td>a<sub>21</sub></td>
 * <td>a<sub>22</sub></td>
 * <td>a<sub>23</sub></td>
 * <td>a<sub>24</sub></td>
 * </tr>
 * <tr>
 * <td>a<sub>31</sub></td>
 * <td>a<sub>32</sub></td>
 * <td>a<sub>33</sub></td>
 * <td>a<sub>34</sub></td>
 * </tr>
 * <tr>
 * <td>a<sub>41</sub></td>
 * <td>a<sub>42</sub></td>
 * <td>a<sub>43</sub></td>
 * <td>a<sub>44</sub></td>
 * </tr>
 * </table>
 * </p>
 * <p>
 * is stored column major in a single array, as follows:
 * </p>
 * <p>
 * <table border="1">
 * <tr>
 * <td>a<sub>11</sub></td>
 * <td>a<sub>21</sub></td>
 * <td>a<sub>31</sub></td>
 * <td>a<sub>41</sub></td>
 * <td>a<sub>12</sub></td>
 * <td>a<sub>22</sub></td>
 * <td>a<sub>32</sub></td>
 * <td>a<sub>42</sub></td>
 * <td>a<sub>13</sub></td>
 * <td>a<sub>23</sub></td>
 * <td>a<sub>33</sub></td>
 * <td>a<sub>43</sub></td>
 * <td>a<sub>14</sub></td>
 * <td>a<sub>24</sub></td>
 * <td>a<sub>34</sub></td>
 * <td>a<sub>44</sub></td>
 * </tr>
 * </table>
 * </p>
 */
public class DenseMatrix extends AbstractDenseMatrix {

    /**
     * Constructor for DenseMatrix
     * 
     * @param r
     *            Reader to get the matrix from
     */
    public DenseMatrix(MatrixVectorReader r) throws IOException {
        // Start with a zero-sized matrix
        super(0, 0);

        // Get matrix information. Use the header if present, else use a safe
        // default
        MatrixInfo info = null;
        if (r.hasInfo())
            info = r.readMatrixInfo();
        else
            info = new MatrixInfo(true, MatrixInfo.MatrixField.Real,
                    MatrixInfo.MatrixSymmetry.General);
        MatrixSize size = r.readMatrixSize(info);

        // Resize the matrix to correct size
        numRows = size.numRows();
        numColumns = size.numColumns();
        data = new double[numRows * numColumns];

        // Check that the matrix is in an acceptable format
        if (info.isPattern())
            throw new UnsupportedOperationException(
                    "Pattern matrices are not supported");
        if (info.isComplex())
            throw new UnsupportedOperationException(
                    "Complex matrices are not supported");

        // Read the entries, in either coordinate or array format
        if (info.isCoordinate()) {

            // Read coordinate data
            int nz = size.numEntries();
            int[] row = new int[nz];
            int[] column = new int[nz];
            double[] entry = new double[nz];
            r.readCoordinate(row, column, entry);

            // Shift indices from 1-offset to 0-offset
            r.add(-1, row);
            r.add(-1, column);

            // Store them
            for (int i = 0; i < nz; ++i)
                set(row[i], column[i], entry[i]);

        } else
            // info.isArray()
            r.readArray(data);

        // Put in missing entries from symmetry or skew symmetry
        if (info.isSymmetric())
            for (int i = 0; i < numRows; ++i)
                for (int j = 0; j < i; ++j)
                    set(j, i, get(i, j));
        else if (info.isSkewSymmetric())
            for (int i = 0; i < numRows; ++i)
                for (int j = 0; j < i; ++j)
                    set(j, i, -get(i, j));
    }

    /**
     * Constructor for DenseMatrix
     * 
     * @param numRows
     *            Number of rows
     * @param numColumns
     *            Number of columns
     */
    public DenseMatrix(int numRows, int numColumns) {
        super(numRows, numColumns);
    }

    /**
     * Constructor for DenseMatrix
     * 
     * @param A
     *            Matrix to copy. A deep copy is made
     */
    public DenseMatrix(Matrix A) {
        super(A);
    }

    /**
     * Constructor for DenseMatrix
     * 
     * @param A
     *            Matrix to copy contents from
     * @param deep
     *            If true, <code>A</code> is copied, else a shallow copy is
     *            made and the matrices share underlying storage. For this,
     *            <code>A</code> must be a dense matrix
     */
    public DenseMatrix(Matrix A, boolean deep) {
        super(A, deep);
    }

    /**
     * Constructor for DenseMatrix. Builds the matrix from a vector
     * 
     * @param x
     *            Vector to copy from. This will form this matrix' single column
     * @param deep
     *            If true, x is copied, if false, the internal storage of this
     *            matrix is the same as that of the vector. In that case,
     *            <code>x</code> must be a <code>DenseVector</code>
     */
    public DenseMatrix(Vector x, boolean deep) {
        super(x.size(), 1);

        if (deep)
            for (VectorEntry e : x)
                set(e.index(), 0, e.get());
        else {
            if (!(x instanceof DenseVector))
                throw new IllegalArgumentException("x must be a DenseVector");
            data = ((DenseVector) x).getData();
        }
    }

    /**
     * Constructor for DenseMatrix. Builds the matrix from a vector
     * 
     * @param x
     *            The vector which forms this matrix' single column. It is
     *            copied, not referenced
     */
    public DenseMatrix(Vector x) {
        this(x, true);
    }

    /**
     * Constructor for DenseMatrix. Builds the matrix from vectors. Each vector
     * will correspond to a column of the matrix
     * 
     * @param x
     *            Vectors which forms the columns of this matrix. Every vector
     *            must have the same size
     */
    public DenseMatrix(Vector[] x) {
        super(x[0].size(), x.length);

        // Ensure correct sizes
        for (Vector v : x)
            if (v.size() != numRows)
                throw new IllegalArgumentException(
                        "All vectors must be of the same size");

        // Copy the contents
        for (int j = 0; j < x.length; ++j)
            for (VectorEntry e : x[j])
                set(e.index(), j, e.get());
    }

    /**
     * Constructor for DenseMatrix. Copies from the passed array
     * 
     * @param values
     *            Arrays to copy from. Every sub-array must have the same size
     */
    public DenseMatrix(double[][] values) {
        super(values.length, values[0].length);

        // Copy the contents
        for (int i = 0; i < values.length; ++i) {
            if (values[i].length != numColumns)
                throw new IllegalArgumentException("Array cannot be jagged");
            for (int j = 0; j < values[i].length; ++j)
                set(i, j, values[i][j]);
        }
    }

    @Override
    public DenseMatrix copy() {
        return new DenseMatrix(this);
    }

    @Override
    void copy(Matrix A) {
        for (MatrixEntry e : A)
            set(e.row(), e.column(), e.get());
    }

    @Override
    public Matrix multAdd(double alpha, Matrix B, Matrix C) {
        if (!(B instanceof DenseMatrix) || !(C instanceof DenseMatrix))
            return super.multAdd(alpha, B, C);

        checkMultAdd(B, C);

        double[] Bd = ((DenseMatrix) B).getData(), Cd = ((DenseMatrix) C)
                .getData();

        BLAS.getInstance().dgemm(Transpose.NoTranspose.netlib(), Transpose.NoTranspose.netlib(),
                C.numRows(), C.numColumns(), numColumns, alpha, data,
                Math.max(1, numRows), Bd, Math.max(1, B.numRows()), 1, Cd,
                Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public Matrix transAmultAdd(double alpha, Matrix B, Matrix C) {
        if (!(B instanceof DenseMatrix) || !(C instanceof DenseMatrix))
            return super.transAmultAdd(alpha, B, C);

        checkTransAmultAdd(B, C);

        double[] Bd = ((DenseMatrix) B).getData(), Cd = ((DenseMatrix) C)
                .getData();

        BLAS.getInstance().dgemm(Transpose.Transpose.netlib(), Transpose.NoTranspose.netlib(),
                C.numRows(), C.numColumns(), numRows, alpha, data,
                Math.max(1, numRows), Bd, Math.max(1, B.numRows()), 1, Cd,
                Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public Matrix transBmultAdd(double alpha, Matrix B, Matrix C) {
        if (!(B instanceof DenseMatrix) || !(C instanceof DenseMatrix))
            return super.transBmultAdd(alpha, B, C);

        checkTransBmultAdd(B, C);

        double[] Bd = ((DenseMatrix) B).getData(), Cd = ((DenseMatrix) C)
                .getData();

        BLAS.getInstance().dgemm(Transpose.NoTranspose.netlib(), Transpose.Transpose.netlib(),
                C.numRows(), C.numColumns(), numColumns, alpha, data,
                Math.max(1, numRows), Bd, Math.max(1, B.numRows()), 1, Cd,
                Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public Matrix transABmultAdd(double alpha, Matrix B, Matrix C) {
        if (!(B instanceof DenseMatrix) || !(C instanceof DenseMatrix))
            return super.transABmultAdd(alpha, B, C);

        checkTransABmultAdd(B, C);

        double[] Bd = ((DenseMatrix) B).getData(), Cd = ((DenseMatrix) C)
                .getData();

        BLAS.getInstance().dgemm(Transpose.Transpose.netlib(), Transpose.Transpose.netlib(),
                C.numRows(), C.numColumns(), numRows, alpha, data,
                Math.max(1, numRows), Bd, Math.max(1, B.numRows()), 1, Cd,
                Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public Matrix rank1(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.rank1(alpha, x, y);

        checkRank1(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dger(numRows, numColumns, alpha, xd, 1, yd, 1, data,
                Math.max(1, numRows));

        return this;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.multAdd(alpha, x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dgemv(Transpose.NoTranspose.netlib(), numRows, numColumns,
                alpha, data, Math.max(numRows, 1), xd, 1, 1, yd, 1);

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.transMultAdd(alpha, x, y);

        checkTransMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        BLAS.getInstance().dgemv(Transpose.Transpose.netlib(), numRows, numColumns, alpha,
                data, Math.max(numRows, 1), xd, 1, 1, yd, 1);

        return y;
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        // We allow non-square matrices, as we then use a least-squares solver
        if (numRows != B.numRows())
            throw new IllegalArgumentException("numRows != B.numRows() ("
                    + numRows + " != " + B.numRows() + ")");
        if (numColumns != X.numRows())
            throw new IllegalArgumentException("numColumns != X.numRows() ("
                    + numColumns + " != " + X.numRows() + ")");
        if (X.numColumns() != B.numColumns())
            throw new IllegalArgumentException(
                    "X.numColumns() != B.numColumns() (" + X.numColumns()
                            + " != " + B.numColumns() + ")");

        if (isSquare())
            return LUsolve(B, X);
        else
            return QRsolve(B, X, Transpose.NoTranspose);
    }

    @Override
    public Vector solve(Vector b, Vector x) {
        DenseMatrix B = new DenseMatrix(b, false), X = new DenseMatrix(x, false);
        solve(B, X);
        return x;
    }

    @Override
    public Matrix transSolve(Matrix B, Matrix X) {
        // We allow non-square matrices, as we then use a least-squares solver
        if (numColumns != B.numRows())
            throw new IllegalArgumentException("numColumns != B.numRows() ("
                    + numColumns + " != " + B.numRows() + ")");
        if (numRows != X.numRows())
            throw new IllegalArgumentException("numRows != X.numRows() ("
                    + numRows + " != " + X.numRows() + ")");
        if (X.numColumns() != B.numColumns())
            throw new IllegalArgumentException(
                    "X.numColumns() != B.numColumns() (" + X.numColumns()
                            + " != " + B.numColumns() + ")");

        return QRsolve(B, X, Transpose.Transpose);
    }

    @Override
    public Vector transSolve(Vector b, Vector x) {
        DenseMatrix B = new DenseMatrix(b, false), X = new DenseMatrix(x, false);
        transSolve(B, X);
        return x;
    }

    Matrix LUsolve(Matrix B, Matrix X) {
        if (!(X instanceof DenseMatrix))
            throw new UnsupportedOperationException("X must be a DenseMatrix");

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        int[] piv = new int[numRows];

        intW info = new intW(0);
        LAPACK.getInstance().dgesv(numRows, B.numColumns(),
                data.clone(), Matrices.ld(numRows), piv, Xd, Matrices.ld(numRows), info);

        if (info.val > 0)
            throw new MatrixSingularException();
        else if (info.val < 0)
            throw new IllegalArgumentException();

        return X;
    }

    Matrix QRsolve(Matrix B, Matrix X, Transpose trans) {
        int nrhs = B.numColumns();

        // Allocate temporary solution matrix
        DenseMatrix Xtmp = new DenseMatrix(Math.max(numRows, numColumns), nrhs);
        int M = trans == Transpose.NoTranspose ? numRows : numColumns;
        for (int j = 0; j < nrhs; ++j)
            for (int i = 0; i < M; ++i)
                Xtmp.set(i, j, B.get(i, j));
        double[] newData = data.clone();

        // Query optimal workspace
        double[] work = new double[1];
        intW info = new intW(0);
        LAPACK.getInstance().dgels(trans.netlib(), numRows, numColumns, nrhs,
                newData, Matrices.ld(numRows), Xtmp.getData(), Matrices.ld(numRows, numColumns),
                work, -1, info);

        // Allocate workspace
        int lwork = -1;
        if (info.val != 0)
            lwork = Math.max(1, Math.min(numRows, numColumns)
                    + Math.max(Math.min(numRows, numColumns), nrhs));
        else
            lwork = Math.max((int) work[0], 1);
        work = new double[lwork];

        // Compute the factorization
        info.val = 0;
        LAPACK.getInstance().dgels(trans.netlib(), numRows, numColumns, nrhs,
                newData, Matrices.ld(numRows), Xtmp.getData(), Matrices.ld(numRows, numColumns),
                work, lwork, info);

        if (info.val < 0)
            throw new IllegalArgumentException();

        // Extract the solution
        int N = trans == Transpose.NoTranspose ? numColumns : numRows;
        for (int j = 0; j < nrhs; ++j)
            for (int i = 0; i < N; ++i)
                X.set(i, j, Xtmp.get(i, j));
        return X;
    }

}
