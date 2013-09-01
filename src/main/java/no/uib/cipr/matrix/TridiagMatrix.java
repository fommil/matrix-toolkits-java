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
import java.util.Iterator;

import com.github.fommil.netlib.LAPACK;
import org.netlib.util.intW;

/**
 * Tridiagonal matrix. Stored in three arrays, one of length <code>n</code>
 * for the diagonal, two of length <code>n-1</code> for the superdiagonal and
 * subdiagonal entries.
 */
public class TridiagMatrix extends AbstractMatrix {

    /**
     * Diagonal, super-diagonal and sub-diagonal
     */
    double[] diag, superDiag, subDiag;

    /**
     * Size of the matrix
     */
    private int n;

    /**
     * Constructor for TridiagMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public TridiagMatrix(int n) {
        super(n, n);

        if (n < 1)
            throw new IllegalArgumentException("n must be >= 1");

        this.n = n;
        diag = new double[n];
        superDiag = new double[n - 1];
        subDiag = new double[n - 1];
    }

    /**
     * Constructor for TridiagMatrix
     * 
     * @param A
     *            Matrix to copy from. Only the tridiagonal part is copied
     */
    public TridiagMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for TridiagMatrix
     * 
     * @param A
     *            Matrix to copy from. Only the tridiagonal part is copied
     * @param deep
     *            True for a deep copy, else it's shallow. For shallow copies,
     *            <code>A</code> must be a <code>TridiagMatrix</code>
     */
    public TridiagMatrix(Matrix A, boolean deep) {
        super(A);

        if (!isSquare())
            throw new IllegalArgumentException(
                    "Tridiagonal matrix must be square");
        if (A.numRows() < 1)
            throw new IllegalArgumentException("numRows must be >= 1");
        n = numRows;

        if (deep) {
            diag = new double[n];
            superDiag = new double[n - 1];
            subDiag = new double[n - 1];
            for (MatrixEntry e : A)
                if (e.row() == e.column() || e.row() == e.column() - 1
                        || e.row() == e.column() + 1)
                    set(e.row(), e.column(), e.get());
        } else {
            TridiagMatrix B = (TridiagMatrix) A;
            this.diag = B.getDiagonal();
            this.subDiag = B.getSubDiagonal();
            this.superDiag = B.getSuperDiagonal();
        }
    }

    /**
     * Returns the diagonal entries. Length equal <code>n</code>
     */
    public double[] getDiagonal() {
        return diag;
    }

    /**
     * Returns the sub diagonal entries. Length equal <code>n-1</code>
     */
    public double[] getSubDiagonal() {
        return subDiag;
    }

    /**
     * Returns the super diagonal entries. Length equal <code>n-1</code>
     */
    public double[] getSuperDiagonal() {
        return superDiag;
    }

    @Override
    public void add(int row, int column, double value) {
        check(row, column);
        if (row == column)
            diag[row] += value;
        else if (row == column + 1)
            subDiag[column] += value;
        else if (row == column - 1)
            superDiag[row] += value;
        else
            throw new IndexOutOfBoundsException(
                    "Insertion index outside of band");
    }

    @Override
    public double get(int row, int column) {
        check(row, column);
        if (row == column)
            return diag[row];
        else if (row == column + 1)
            return subDiag[column];
        else if (row == column - 1)
            return superDiag[row];
        else
            return 0;
    }

    @Override
    public void set(int row, int column, double value) {
        check(row, column);
        if (row == column)
            diag[row] = value;
        else if (row == column + 1)
            subDiag[column] = value;
        else if (row == column - 1)
            superDiag[row] = value;
        else
            throw new IndexOutOfBoundsException(
                    "Insertion index outside of band");
    }

    @Override
    public TridiagMatrix copy() {
        return new TridiagMatrix(this);
    }

    @Override
    public TridiagMatrix zero() {
        Arrays.fill(diag, 0);
        Arrays.fill(subDiag, 0);
        Arrays.fill(superDiag, 0);
        return this;
    }

    @Override
    public Matrix solve(Matrix B, Matrix X) {
        if (!(X instanceof DenseMatrix))
            throw new UnsupportedOperationException("X must be a DenseMatrix");

        checkSolve(B, X);

        double[] Xd = ((DenseMatrix) X).getData();

        X.set(B);

        intW info = new intW(0);
        LAPACK.getInstance().dgtsv(numRows, X.numColumns(),
                subDiag.clone(), diag.clone(), superDiag.clone(), Xd, Matrices.ld(numRows), info);

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
        double[] otherDiag = subDiag;
        subDiag = superDiag;
        superDiag = otherDiag;
        return this;
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new TridiagMatrixIterator();
    }

    /**
     * Iterator over a tridiagonal matrix
     */
    private class TridiagMatrixIterator extends RefMatrixIterator {

        /**
         * Current band, starting with the main diagonal
         */
        private double[] band = diag;

        /**
         * Band index
         */
        private int bandIndex;

        /**
         * Which band in use (0 for main, 1 for sub, 2 for super)
         */
        private int whichBand;

        @Override
        public boolean hasNext() {
            return whichBand < 3;
        }

        @Override
        public MatrixEntry next() {
            entry.update(row, column);

            // Move in the band
            if (bandIndex < band.length - 1)
                bandIndex++;
            else {
                // Move to the next band
                bandIndex = 0;
                whichBand++;

                if (whichBand == 1)
                    band = subDiag;
                else if (whichBand == 2)
                    band = superDiag;

                // If the off-diagonals are zero-sized, we are done
                // This happens if the matrix is 1*1
                if (band.length == 0)
                    whichBand = 3;
            }

            // Set row index
            if (whichBand == 1)
                row = bandIndex + 1;
            else
                row = bandIndex;

            // Set column index
            if (whichBand == 2)
                column = bandIndex + 1;
            else
                column = bandIndex;

            return entry;
        }

    }

}
