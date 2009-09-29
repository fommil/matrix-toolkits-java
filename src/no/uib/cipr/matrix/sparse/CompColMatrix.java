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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.io.MatrixInfo;
import no.uib.cipr.matrix.io.MatrixSize;
import no.uib.cipr.matrix.io.MatrixVectorReader;

/**
 * Compressed column storage (CCS) matrix
 */
public class CompColMatrix extends AbstractMatrix {

    /**
     * Matrix data
     */
    double[] data;

    /**
     * Column indices. These are kept sorted within each row.
     */
    int[] columnPointer;

    /**
     * Indices to the start of each row
     */
    int[] rowIndex;

    /**
     * Constructor for CompColMatrix
     * 
     * @param r
     *            Reader to get sparse matrix from
     */
    public CompColMatrix(MatrixVectorReader r) throws IOException {
        // Start with a zero-sized matrix
        super(0, 0);

        // Get matrix information. Use the header if present, else just assume
        // that the matrix stores real numbers without any symmetry
        MatrixInfo info = null;
        if (r.hasInfo())
            info = r.readMatrixInfo();
        else
            info = new MatrixInfo(true, MatrixInfo.MatrixField.Real,
                    MatrixInfo.MatrixSymmetry.General);

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

        // Resize the matrix to correct size
        MatrixSize size = r.readMatrixSize(info);
        numRows = size.numRows();
        numColumns = size.numColumns();

        // Start reading entries
        int numEntries = size.numEntries();
        int[] row = new int[numEntries];
        int[] column = new int[numEntries];
        double[] entry = new double[numEntries];
        r.readCoordinate(row, column, entry);

        // Shift the indices from 1 based to 0 based
        r.add(-1, row);
        r.add(-1, column);

        // Find the number of entries on each column
        List<Set<Integer>> cnz = new ArrayList<Set<Integer>>(numColumns);
        for (int i = 0; i < numColumns; ++i)
            cnz.add(new HashSet<Integer>());

        for (int i = 0; i < numEntries; ++i)
            cnz.get(column[i]).add(row[i]);

        // Allocate some more in case of symmetry
        if (info.isSymmetric() || info.isSkewSymmetric())
            for (int i = 0; i < numEntries; ++i)
                if (row[i] != column[i])
                    cnz.get(row[i]).add(column[i]);

        int[][] nz = new int[numColumns][];
        for (int i = 0; i < numColumns; ++i) {
            nz[i] = new int[cnz.get(i).size()];
            int j = 0;
            for (Integer rowind : cnz.get(i))
                nz[i][j++] = rowind;
        }

        // Create the sparse matrix structure
        construct(nz);

        // Insert the entries
        for (int i = 0; i < size.numEntries(); ++i)
            set(row[i], column[i], entry[i]);

        // Put in extra entries from symmetry or skew symmetry
        if (info.isSymmetric())
            for (int i = 0; i < numEntries; ++i) {
                if (row[i] != column[i])
                    set(column[i], row[i], entry[i]);
            }
        else if (info.isSkewSymmetric())
            for (int i = 0; i < numEntries; ++i) {
                if (row[i] != column[i])
                    set(column[i], row[i], -entry[i]);
            }
    }

    /**
     * Constructor for CompColMatrix
     * 
     * @param numRows
     *            Number of rows
     * @param numColumns
     *            Number of columns
     * @param nz
     *            The nonzero column indices on each column
     */
    public CompColMatrix(int numRows, int numColumns, int[][] nz) {
        super(numRows, numColumns);
        construct(nz);
    }

    private void construct(int[][] nz) {
        int nnz = 0;
        for (int i = 0; i < nz.length; ++i)
            nnz += nz[i].length;

        columnPointer = new int[numColumns + 1];
        rowIndex = new int[nnz];
        data = new double[nnz];

        if (nz.length != numColumns)
            throw new IllegalArgumentException("nz.length != numColumns");

        for (int i = 1; i <= numColumns; ++i) {
            columnPointer[i] = columnPointer[i - 1] + nz[i - 1].length;

            for (int j = columnPointer[i - 1], k = 0; j < columnPointer[i]; ++j, ++k) {
                rowIndex[j] = nz[i - 1][k];
                if (nz[i - 1][k] < 0 || nz[i - 1][k] >= numRows)
                    throw new IllegalArgumentException("nz[" + (i - 1) + "]["
                            + k + "]=" + nz[i - 1][k]
                            + ", which is not a valid row index");
            }

            Arrays.sort(rowIndex, columnPointer[i - 1], columnPointer[i]);
        }
    }

    private void construct(Matrix A, boolean deep) {
        if (deep) {
            if (A instanceof CompColMatrix) {
                CompColMatrix Ac = (CompColMatrix) A;
                data = new double[Ac.data.length];
                columnPointer = new int[Ac.columnPointer.length];
                rowIndex = new int[Ac.rowIndex.length];

                System.arraycopy(Ac.data, 0, data, 0, data.length);
                System.arraycopy(Ac.columnPointer, 0, columnPointer, 0,
                        columnPointer.length);
                System.arraycopy(Ac.rowIndex, 0, rowIndex, 0, rowIndex.length);
            } else {

                List<Set<Integer>> cnz = new ArrayList<Set<Integer>>(numColumns);
                for (int i = 0; i < numColumns; ++i)
                    cnz.add(new HashSet<Integer>());

                for (MatrixEntry e : A)
                    cnz.get(e.column()).add(e.row());

                int[][] nz = new int[numColumns][];
                for (int i = 0; i < numColumns; ++i) {
                    nz[i] = new int[cnz.get(i).size()];
                    int j = 0;
                    for (Integer rowind : cnz.get(i))
                        nz[i][j++] = rowind;
                }

                construct(nz);
                set(A);

            }
        } else {
            CompColMatrix Ac = (CompColMatrix) A;
            columnPointer = Ac.getColumnPointers();
            rowIndex = Ac.getRowIndices();
            data = Ac.getData();
        }
    }

    /**
     * Constructor for CompColMatrix
     * 
     * @param A
     *            Copies from this matrix
     * @param deep
     *            True if the copy is to be deep. If it is a shallow copy,
     *            <code>A</code> must be a <code>CompColMatrix</code>
     */
    public CompColMatrix(Matrix A, boolean deep) {
        super(A);
        construct(A, deep);
    }

    /**
     * Constructor for CompColMatrix
     * 
     * @param A
     *            Copies from this matrix. The copy will be deep
     */
    public CompColMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Returns the column pointers
     */
    public int[] getColumnPointers() {
        return columnPointer;
    }

    /**
     * Returns the row indices
     */
    public int[] getRowIndices() {
        return rowIndex;
    }

    /**
     * Returns the internal data storage
     */
    public double[] getData() {
        return data;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.multAdd(alpha, x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        // y = 1/alpha * y
        y.scale(1 / alpha);

        // y = A*x + y
        for (int i = 0; i < numColumns; ++i)
            for (int j = columnPointer[i]; j < columnPointer[i + 1]; ++j)
                yd[rowIndex[j]] += data[j] * xd[i];

        // y = alpha*y = alpha*A*x + y
        return y.scale(alpha);
    }

    @Override
    public Vector transMult(Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.transMult(x, y);

        checkTransMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData();
        double[] yd = ((DenseVector) y).getData();

        for (int i = 0; i < numColumns; ++i) {
            double dot = 0;
            for (int j = columnPointer[i]; j < columnPointer[i + 1]; ++j)
                dot += data[j] * xd[rowIndex[j]];
            yd[i] = dot;
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

        for (int i = 0; i < numColumns; ++i) {
            double dot = 0;
            for (int j = columnPointer[i]; j < columnPointer[i + 1]; ++j)
                dot += data[j] * xd[rowIndex[j]];
            yd[i] += alpha * dot;
        }

        return y;
    }

    @Override
    public void set(int row, int column, double value) {
        check(row, column);

        int index = getIndex(row, column);
        data[index] = value;
    }

    @Override
    public void add(int row, int column, double value) {
        check(row, column);

        int index = getIndex(row, column);
        data[index] += value;
    }

    @Override
    public double get(int row, int column) {
        check(row, column);

        int index = no.uib.cipr.matrix.sparse.Arrays.binarySearch(rowIndex,
                row, columnPointer[column], columnPointer[column + 1]);

        if (index >= 0)
            return data[index];
        else
            return 0;
    }

    /**
     * Finds the insertion index
     */
    private int getIndex(int row, int column) {
        int i = no.uib.cipr.matrix.sparse.Arrays.binarySearch(rowIndex, row,
                columnPointer[column], columnPointer[column + 1]);

        if (i != -1 && rowIndex[i] == row)
            return i;
        else
            throw new IndexOutOfBoundsException("Entry (" + (row + 1) + ", "
                    + (column + 1) + ") is not in the matrix structure");
    }

    @Override
    public CompColMatrix copy() {
        return new CompColMatrix(this);
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new CompColMatrixIterator();
    }

    @Override
    public CompColMatrix zero() {
        Arrays.fill(data, 0);
        return this;
    }

    /**
     * Iterator over a compressed column matrix
     */
    private class CompColMatrixIterator implements Iterator<MatrixEntry> {

        private int column, cursor;

        private CompColMatrixEntry entry = new CompColMatrixEntry();

        public CompColMatrixIterator() {
            // Find first non-empty column
            nextNonEmptyColumn();
        }

        /**
         * Locates the first non-empty column, starting at the current. After
         * the new column has been found, the cursor is also updated
         */
        private void nextNonEmptyColumn() {
            while (column < numColumns()
                    && columnPointer[column] == columnPointer[column + 1])
                column++;
            cursor = columnPointer[column];
        }

        public boolean hasNext() {
            return cursor < data.length;
        }

        public MatrixEntry next() {
            entry.update(column, cursor);

            // Next position is in the same column
            if (cursor < columnPointer[column + 1] - 1)
                cursor++;

            // Next position is at the following (non-empty) column
            else {
                column++;
                nextNonEmptyColumn();
            }

            return entry;
        }

        public void remove() {
            entry.set(0);
        }

    }

    /**
     * Entry of a compressed column matrix
     */
    private class CompColMatrixEntry implements MatrixEntry {

        private int column, cursor;

        /**
         * Updates the entry
         */
        public void update(int column, int cursor) {
            this.column = column;
            this.cursor = cursor;
        }

        public int row() {
            return rowIndex[cursor];
        }

        public int column() {
            return column;
        }

        public double get() {
            return data[cursor];
        }

        public void set(double value) {
            data[cursor] = value;
        }
    }
}
