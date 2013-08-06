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

import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.io.MatrixInfo;
import no.uib.cipr.matrix.io.MatrixSize;
import no.uib.cipr.matrix.io.MatrixVectorReader;

/**
 * Compressed diagonal storage (CDS) matrix
 */
public class CompDiagMatrix extends AbstractMatrix {

    /**
     * The diagonals
     */
    double[][] diag;

    /**
     * Indices to the start of the diagonal, relative to the main diagonal.
     * Positive means the number of diagonals shifted up, while negative is the
     * number of diagonals shifted down
     */
    int[] ind;

    /**
     * Constructor for CompDiagMatrix
     * 
     * @param r
     *            Reader to get sparse matrix from
     */
    public CompDiagMatrix(MatrixVectorReader r) throws IOException {
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

        // Check that the matrix is in an acceptable format
        if (info.isPattern())
            throw new UnsupportedOperationException(
                    "Pattern matrices are not supported");
        if (info.isDense())
            throw new UnsupportedOperationException(
                    "Dense matrices are not supported");
        if (info.isComplex())
            throw new UnsupportedOperationException(
                    "Complex matrices are not supported");

        // Start reading entries
        int[] row = new int[size.numEntries()], column = new int[size
                .numEntries()];
        double[] entry = new double[size.numEntries()];
        r.readCoordinate(row, column, entry);

        // Shift the indices from 1 based to 0 based
        r.add(-1, row);
        r.add(-1, column);

        // Find all the diagonals so that we can preallocate
        Set<Integer> diags = new TreeSet<Integer>();
        for (int i = 0; i < size.numEntries(); ++i)
            diags.add(getDiagonal(row[i], column[i]));

        if (info.isSymmetric() || info.isSkewSymmetric())
            for (int i = 0; i < size.numEntries(); ++i)
                if (row[i] != column[i])
                    diags.add(getDiagonal(column[i], row[i]));

        // Convert into an integer array
        int[] ind = new int[diags.size()];
        {
            Integer[] ints = new Integer[diags.size()];
            diags.toArray(ints);
            for (int i = 0; i < diags.size(); ++i)
                ind[i] = ints[i];
        }

        // Create the structure with preallocation
        construct(ind);

        // Insert the entries
        for (int i = 0; i < size.numEntries(); ++i)
            set(row[i], column[i], entry[i]);

        // Put in missing entries from symmetry or skew symmetry
        if (info.isSymmetric())
            for (int i = 0; i < size.numEntries(); ++i) {
                if (row[i] != column[i])
                    set(column[i], row[i], entry[i]);
            }
        else if (info.isSkewSymmetric())
            for (int i = 0; i < size.numEntries(); ++i) {
                if (row[i] != column[i])
                    set(column[i], row[i], -entry[i]);
            }
    }

    /**
     * Creates a new sparse matrix with the given diagonals preallocated. A
     * negative index is a subdiagonal, positive is superdiagonal
     */
    public CompDiagMatrix(int numRows, int numColumns, int[] diagonal) {
        super(numRows, numColumns);
        construct(diagonal);
    }

    private void construct(int[] diagonal) {
        diag = new double[diagonal.length][];
        ind = new int[diagonal.length];

        // Keep the diagonal indices sorted
        int[] sortedDiagonal = new int[diagonal.length];
        System.arraycopy(diagonal, 0, sortedDiagonal, 0, diagonal.length);
        Arrays.sort(sortedDiagonal);

        for (int i = 0; i < diagonal.length; ++i) {
            ind[i] = sortedDiagonal[i];
            diag[i] = new double[getDiagSize(sortedDiagonal[i])];
        }
    }

    /**
     * Creates a new sparse matrix without preallocation
     */
    public CompDiagMatrix(int numRows, int numColumns) {
        this(numRows, numColumns, new int[0]);
    }

    /**
     * Creates a new sparse matrix copied from the given matrix. Can take a deep
     * copy or a shallow copy. For the latter, the supplied matrix must be a
     * CompDiagMatrix. Preallocation is also possible, but is only used for the
     * deep copy.
     */
    public CompDiagMatrix(Matrix A, int[] diagonal, boolean deep) {
        super(A);

        if (deep) {
            construct(diagonal);
            set(A);
        } else {
            CompDiagMatrix Ac = (CompDiagMatrix) A;
            diag = Ac.getDiagonals();
            ind = Ac.getIndex();
        }
    }

    /**
     * Creates a new sparse matrix copied from the given matrix. Takes a deep
     * copy, with possibility to specify preallocation
     */
    public CompDiagMatrix(Matrix A, int[] diagonal) {
        this(A, diagonal, true);
    }

    /**
     * Creates a new sparse matrix copied from the given matrix. Can take a deep
     * copy or a shallow copy. For the latter, the supplied matrix must be a
     * CompDiagMatrix. No preallocation is done
     */
    public CompDiagMatrix(Matrix A, boolean deep) {
        this(A, new int[0], deep);
    }

    /**
     * Creates a new sparse matrix copied from the given matrix. Takes a deep
     * copy without preallocation
     */
    public CompDiagMatrix(Matrix A) {
        this(A, new int[0], true);
    }

    /**
     * Returns the internal diagonal storage
     */
    public double[][] getDiagonals() {
        return diag;
    }

    /**
     * Returns the diagonal offsets
     */
    public int[] getIndex() {
        return ind;
    }

    @Override
    public void add(int row, int column, double value) {
        check(row, column);

        int diagonal = getCompDiagIndex(row, column);

        diag[diagonal][getOnDiagIndex(row, column)] += value;
    }

    @Override
    public double get(int row, int column) {
        check(row, column);

        int diagonal = Arrays.binarySearch(ind, getDiagonal(row, column));

        if (diagonal >= 0)
            return diag[diagonal][getOnDiagIndex(row, column)];
        else
            return 0;
    }

    @Override
    public void set(int row, int column, double value) {
        check(row, column);

        int diagonal = getCompDiagIndex(row, column);

        diag[diagonal][getOnDiagIndex(row, column)] = value;
    }

    private int getDiagonal(int row, int column) {
        return column - row;
    }

    private int getOnDiagIndex(int row, int column) {
        return row > column ? column : row;
    }

    private int getCompDiagIndex(int row, int column) {
        int diagonal = getDiagonal(row, column);

        // Check if the diagonal is already present
        int index = no.uib.cipr.matrix.sparse.Arrays.binarySearchGreater(ind,
                diagonal);
        if (index < ind.length && ind[index] == diagonal)
            return index;

        // Need to allocate new diagonal. Get the diagonal size
        int size = getDiagSize(diagonal);

        // Allocate new primary structure
        double[] newDiag = new double[size];
        double[][] newDiagArray = new double[diag.length + 1][];
        int[] newInd = new int[ind.length + 1];

        // Move data from the old into the new structure
        System.arraycopy(ind, 0, newInd, 0, index);
        System.arraycopy(ind, index, newInd, index + 1, ind.length - index);
        for (int i = 0; i < index; ++i)
            newDiagArray[i] = diag[i];
        for (int i = index; i < diag.length; ++i)
            newDiagArray[i + 1] = diag[i];

        newInd[index] = diagonal;
        newDiagArray[index] = newDiag;

        // Update pointers
        ind = newInd;
        diag = newDiagArray;

        return index;
    }

    /**
     * Finds the size of the requested diagonal to be allocated
     */
    private int getDiagSize(int diagonal) {
        if (diagonal < 0)
            return Math.min(numRows + diagonal, numColumns);
        else
            return Math.min(numRows, numColumns - diagonal);
    }

    @Override
    public Matrix copy() {
        return new CompDiagMatrix(this, ind);
    }

    @Override
    public Matrix zero() {
        for (int i = 0; i < diag.length; ++i)
            Arrays.fill(diag[i], 0);
        return this;
    }

    @Override
    public Vector mult(Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.mult(x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData();
        double[] yd = ((DenseVector) y).getData();

        y.zero();

        for (int i = 0; i < ind.length; ++i) {
            int row = ind[i] < 0 ? -ind[i] : 0;
            int column = ind[i] > 0 ? ind[i] : 0;
            double[] locDiag = diag[i];
            for (int j = 0; j < locDiag.length; ++j, ++row, ++column)
                yd[row] += locDiag[j] * xd[column];
        }

        return y;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.multAdd(alpha, x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData();
        double[] yd = ((DenseVector) y).getData();

        for (int i = 0; i < ind.length; ++i) {
            int row = ind[i] < 0 ? -ind[i] : 0;
            int column = ind[i] > 0 ? ind[i] : 0;
            double[] locDiag = diag[i];
            for (int j = 0; j < locDiag.length; ++j, ++row, ++column)
                yd[row] += alpha * locDiag[j] * xd[column];
        }

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.transMultAdd(alpha, x, y);

        checkTransMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData();
        double[] yd = ((DenseVector) y).getData();

        for (int i = 0; i < ind.length; ++i) {
            int row = ind[i] < 0 ? -ind[i] : 0;
            int column = ind[i] > 0 ? ind[i] : 0;
            double[] locDiag = diag[i];
            for (int j = 0; j < locDiag.length; ++j, ++row, ++column)
                yd[column] += alpha * locDiag[j] * xd[row];
        }

        return y;
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new CompDiagMatrixIterator();
    }

    /**
     * Iterator over a compressed diagonal matrix
     */
    private class CompDiagMatrixIterator implements Iterator<MatrixEntry> {

        private int diagonal, index;

        private CompDiagMatrixEntry entry = new CompDiagMatrixEntry();

        public boolean hasNext() {
            return diagonal < diag.length;
        }

        public MatrixEntry next() {
            entry.update(diagonal, index);

            // Move along current diagonal
            if (index < diag[diagonal].length - 1)
                index++;

            // Move to the next diagonal
            else {
                diagonal++;
                index = 0;
            }

            return entry;
        }

        public void remove() {
            entry.set(0);
        }

    }

    /**
     * Entry of a compressed diagonal matrix
     */
    private class CompDiagMatrixEntry implements MatrixEntry {

        private int diagonal, index;

        public void update(int diagonal, int index) {
            this.diagonal = diagonal;
            this.index = index;
        }

        public int row() {
            return index + (ind[diagonal] < 0 ? -ind[diagonal] : 0);
        }

        public int column() {
            return index + (ind[diagonal] > 0 ? ind[diagonal] : 0);
        }

        public double get() {
            return diag[diagonal][index];
        }

        public void set(double value) {
            diag[diagonal][index] = value;
        }

    }

}
