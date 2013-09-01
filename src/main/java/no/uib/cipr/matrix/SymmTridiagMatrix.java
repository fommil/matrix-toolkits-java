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
 * Symmetrical tridiagonal matrix. Storage as for
 * {@link no.uib.cipr.matrix.TridiagMatrix TridiagMatrix}, but only one
 * off-diagonal array is stored.
 */
public class SymmTridiagMatrix extends AbstractMatrix {

    /**
     * Diagonal and off-diagonal
     */
    double[] diag, offDiag;

    /**
     * Size of the matrix
     */
    int n;

    /**
     * Constructor for SymmTridiagMatrix
     * 
     * @param diag
     *            Main diagonal
     * @param offDiag
     *            Offdiagonals, both upper and lower
     * @param n
     *            Size of the matrix. The main diagonal must be at least as long
     *            as n, and the off diagonal part must be at least as long as
     *            n-1
     */
    public SymmTridiagMatrix(double[] diag, double[] offDiag, int n) {
        super(n, n);

        this.n = n;
        if (n < 1)
            throw new IllegalArgumentException("n must be >= 1");

        if (diag.length < n)
            throw new IllegalArgumentException("diag.length < n");
        if (offDiag.length < n - 1)
            throw new IllegalArgumentException("offDiag.length < n - 1");

        this.diag = diag;
        this.offDiag = offDiag;
    }

    /**
     * Constructor for SymmTridiagMatrix
     * 
     * @param diag
     *            Main diagonal
     * @param offDiag
     *            Offdiagonals. Must be one shorter than diag
     */
    public SymmTridiagMatrix(double[] diag, double[] offDiag) {
        this(diag, offDiag, diag.length);
    }

    /**
     * Constructor for SymmTridiagMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns. <code>n</code>
     *            cannot be zero
     */
    public SymmTridiagMatrix(int n) {
        super(n, n);

        if (n < 1)
            throw new IllegalArgumentException("n must be >= 1");

        this.n = numRows;
        diag = new double[n];
        offDiag = new double[n - 1];
    }

    /**
     * Constructor for SymmTridiagMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only main and the superdiagonal
     *            is copied over
     */
    public SymmTridiagMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for SymmTridiagMatrix
     * 
     * @param A
     *            Matrix to copy contents from. Only main and the superdiagonal
     *            is copied over. It must be square and cannot have any zero
     *            dimension lengths
     * @param deep
     *            True for a deep copy. For shallow copies <code>A</code> must
     *            be a <code>SymmTridiagMatrix</code>
     */
    public SymmTridiagMatrix(Matrix A, boolean deep) {
        super(A);

        if (!isSquare())
            throw new IllegalArgumentException(
                    "Symmetric matrix must be square");
        if (A.numRows() < 1)
            throw new IllegalArgumentException("numRows must be >= 1");
        n = numRows;

        if (deep) {
            diag = new double[n];
            offDiag = new double[Math.max(n - 1, 0)];
            for (MatrixEntry e : A)
                if (e.row() == e.column() || e.row() == e.column() + 1)
                    set(e.row(), e.column(), e.get());
        } else {
            SymmTridiagMatrix B = (SymmTridiagMatrix) A;
            this.diag = B.getDiagonal();
            this.offDiag = B.getOffDiagonal();
        }
    }

    /**
     * Returns the diagonal entries. Length equal <code>n</code>
     */
    public double[] getDiagonal() {
        return diag;
    }

    /**
     * Returns the off diagonal entries. Length equal <code>n-1</code>
     */
    public double[] getOffDiagonal() {
        return offDiag;
    }

    @Override
    public void add(int row, int column, double value) {
        check(row, column);
        if (row == column)
            diag[row] += value;
        else if (row == column + 1)
            offDiag[column] += value;
        else if (row != column - 1)
            throw new IndexOutOfBoundsException(
                    "Insertion index outside of band");
    }

    @Override
    public double get(int row, int column) {
        check(row, column);
        if (row == column)
            return diag[row];
        else if (row == column + 1)
            return offDiag[column];
        else if (row == column - 1)
            return offDiag[row];
        else
            return 0;
    }

    @Override
    public void set(int row, int column, double value) {
        check(row, column);
        if (row == column)
            diag[row] = value;
        else if (row == column + 1)
            offDiag[column] = value;
        else if (row != column - 1)
            throw new IndexOutOfBoundsException(
                    "Insertion index outside of band");
    }

    @Override
    public SymmTridiagMatrix copy() {
        return new SymmTridiagMatrix(this);
    }

    @Override
    public SymmTridiagMatrix zero() {
        Arrays.fill(diag, 0);
        Arrays.fill(offDiag, 0);
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
                offDiag.clone(), diag.clone(), offDiag.clone(), Xd, Matrices.ld(numRows), info);

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
    public Matrix transSolve(Matrix B, Matrix X) {
        return solve(B, X);
    }

    @Override
    public Vector transSolve(Vector b, Vector x) {
        return solve(b, x);
    }

    @Override
    public Matrix transpose() {
        return this;
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new SymmTridiagMatrixIterator();
    }

    /**
     * Iterator over a symmetrical, tridiagonal matrix
     */
    private class SymmTridiagMatrixIterator extends RefMatrixIterator {

        /**
         * Current band, starting with the main diagonal
         */
        private double[] band = diag;

        /**
         * Band index
         */
        private int bandIndex;

        /**
         * Which band in use (0 for main, 1 for off)
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
                // Move to the off-diagonal (twice)
                bandIndex = 0;
                whichBand++;

                band = offDiag;

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
