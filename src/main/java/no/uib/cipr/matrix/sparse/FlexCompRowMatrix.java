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

import java.util.Arrays;
import java.util.Iterator;

import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;
import no.uib.cipr.matrix.sparse.SuperIterator.SuperIteratorEntry;

/**
 * Matrix stored row-wise into sparse vectors
 */
public class FlexCompRowMatrix extends AbstractMatrix {

    /**
     * Matrix data
     */
    SparseVector[] rowD;

    /**
     * Constructor for FlexCompRowMatrix
     * 
     * @param numRows
     *            Number of rows
     * @param numColumns
     *            Number of column
     */
    public FlexCompRowMatrix(int numRows, int numColumns) {
        super(numRows, numColumns);

        rowD = new SparseVector[numRows];
        for (int i = 0; i < numRows; ++i)
            rowD[i] = new SparseVector(numColumns);
    }

    /**
     * Constructor for FlexCompRowMatrix
     * 
     * @param A
     *            Matrix to copy contents from
     * @param deep
     *            True for a deep copy, false for a reference copy. A reference
     *            copy can only be made of an <code>FlexCompRowMatrix</code>
     */
    public FlexCompRowMatrix(Matrix A, boolean deep) {
        super(A);
        rowD = new SparseVector[numRows];

        if (deep) {
            for (int i = 0; i < numRows; ++i)
                rowD[i] = new SparseVector(numColumns);
            set(A);
        } else {
            FlexCompRowMatrix Ar = (FlexCompRowMatrix) A;
            for (int i = 0; i < numRows; ++i)
                rowD[i] = Ar.getRow(i);
        }
    }

    /**
     * Constructor for FlexCompRowMatrix
     * 
     * @param A
     *            Matrix to copy contents from. The copy will be deep
     */
    public FlexCompRowMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Returns the given row
     */
    public SparseVector getRow(int i) {
        return rowD[i];
    }

    /**
     * Sets the given row equal the passed vector
     */
    public void setRow(int i, SparseVector x) {
        if (x.size() != numColumns)
            throw new IllegalArgumentException(
                    "New row must be of the same size as existing row");
        rowD[i] = x;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        checkMultAdd(x, y);

        for (int i = 0; i < numRows; ++i)
            y.add(i, alpha * rowD[i].dot(x));

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.transMultAdd(alpha, x, y);

        checkTransMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData(), yd = ((DenseVector) y)
                .getData();

        // y = 1/alpha * y
        y.scale(1. / alpha);

        // y = A'x + y
        for (int i = 0; i < numRows; ++i) {
            SparseVector v = rowD[i];
            int[] index = v.getIndex();
            double[] data = v.getData();
            int length = v.getUsed();
            for (int j = 0; j < length; ++j)
                yd[index[j]] += data[j] * xd[i];
        }

        // y = alpha*y = alpha * A'x + y
        return y.scale(alpha);
    }

    @Override
    public void add(int row, int column, double value) {
        rowD[row].add(column, value);
    }

    @Override
    public void set(int row, int column, double value) {
        rowD[row].set(column, value);
    }

    @Override
    public double get(int row, int column) {
        return rowD[row].get(column);
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new RowMatrixIterator();
    }

    @Override
    public Matrix copy() {
        return new FlexCompRowMatrix(this);
    }

    @Override
    public FlexCompRowMatrix zero() {
        for (int i = 0; i < numRows; ++i)
            rowD[i].zero();
        return this;
    }

    @Override
    public Matrix set(Matrix B) {
        if (!(B instanceof FlexCompRowMatrix))
            return super.set(B);

        checkSize(B);

        FlexCompRowMatrix Bc = (FlexCompRowMatrix) B;

        for (int i = 0; i < numRows; ++i)
            rowD[i].set(Bc.rowD[i]);

        return this;
    }

    /**
     * Tries to store the matrix as compactly as possible
     */
    public void compact() {
        for (Vector v : rowD)
            ((SparseVector) v).compact();
    }

    /**
     * Iterator over a matrix stored vectorwise by rows
     */
    private class RowMatrixIterator implements Iterator<MatrixEntry> {

        /**
         * Iterates over each row vector
         */
        private SuperIterator<SparseVector, VectorEntry> iterator = new SuperIterator<SparseVector, VectorEntry>(
                Arrays.asList(rowD));

        /**
         * Entry returned
         */
        private RowMatrixEntry entry = new RowMatrixEntry();

        public boolean hasNext() {
            return iterator.hasNext();
        }

        public MatrixEntry next() {
            SuperIteratorEntry<VectorEntry> se = iterator.next();
            entry.update(se.index(), se.get());
            return entry;
        }

        public void remove() {
            iterator.remove();
        }

    }

    /**
     * Entry of a matrix stored vectorwise by rows
     */
    private static class RowMatrixEntry implements MatrixEntry {

        private int row;

        private VectorEntry entry;

        public void update(int row, VectorEntry entry) {
            this.row = row;
            this.entry = entry;
        }

        public int row() {
            return row;
        }

        public int column() {
            return entry.index();
        }

        public double get() {
            return entry.get();
        }

        public void set(double value) {
            entry.set(value);
        }

    }

}
