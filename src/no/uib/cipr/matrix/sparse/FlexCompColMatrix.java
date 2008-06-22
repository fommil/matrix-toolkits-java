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
 * Matrix stored column-wise into sparse vectors
 */
public class FlexCompColMatrix extends AbstractMatrix {

    /**
     * Matrix data
     */
    SparseVector[] colD;

    /**
     * Constructor for FlexCompColMatrix
     * 
     * @param numRows
     *            Number of rows
     * @param numColumns
     *            Number of column
     */
    public FlexCompColMatrix(int numRows, int numColumns) {
        super(numRows, numColumns);

        colD = new SparseVector[numColumns];
        for (int i = 0; i < numColumns; ++i)
            colD[i] = new SparseVector(numRows);
    }

    /**
     * Constructor for FlexCompColMatrix
     * 
     * @param A
     *            Matrix to copy contents from
     * @param deep
     *            True for a deep copy, false for a reference copy. A reference
     *            copy can only be made of an <code>FlexCompColMatrix</code>
     */
    public FlexCompColMatrix(Matrix A, boolean deep) {
        super(A);
        colD = new SparseVector[numColumns];

        if (deep) {
            for (int i = 0; i < numColumns; ++i)
                colD[i] = new SparseVector(numRows);
            set(A);
        } else {
            FlexCompColMatrix Ar = (FlexCompColMatrix) A;
            for (int i = 0; i < numColumns; ++i)
                colD[i] = Ar.getColumn(i);
        }
    }

    /**
     * Constructor for FlexCompColMatrix
     * 
     * @param A
     *            Matrix to copy contents from. The copy will be deep
     */
    public FlexCompColMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Returns the given column
     */
    public SparseVector getColumn(int i) {
        return colD[i];
    }

    /**
     * Sets the given column equal the passed vector
     */
    public void setColumn(int i, SparseVector x) {
        if (x.size() != numRows)
            throw new IllegalArgumentException(
                    "New column must be of the same size as existing column");
        colD[i] = x;
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DenseVector) || !(y instanceof DenseVector))
            return super.multAdd(alpha, x, y);

        checkMultAdd(x, y);

        double[] xd = ((DenseVector) x).getData();
        double[] yd = ((DenseVector) y).getData();

        // y = 1/alpha * y
        y.scale(1. / alpha);

        // y = A*x + y
        for (int i = 0; i < numColumns; ++i) {
            SparseVector v = colD[i];
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
    public Vector transMultAdd(final double alpha, final Vector x,
            final Vector y) {
        checkTransMultAdd(x, y);

        for (int i = 0; i < numColumns; ++i)
            y.add(i, alpha * colD[i].dot(x));

        return y;
    }

    @Override
    public void add(int row, int column, double value) {
        colD[column].add(row, value);
    }

    @Override
    public void set(int row, int column, double value) {
        colD[column].set(row, value);
    }

    @Override
    public double get(int row, int column) {
        return colD[column].get(row);
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new ColMatrixIterator();
    }

    @Override
    public FlexCompColMatrix copy() {
        return new FlexCompColMatrix(this);
    }

    @Override
    public FlexCompColMatrix zero() {
        for (int i = 0; i < numColumns; ++i)
            colD[i].zero();
        return this;
    }

    /**
     * Tries to store the matrix as compactly as possible
     */
    public void compact() {
        for (Vector v : colD)
            ((SparseVector) v).compact();
    }

    /**
     * Iterator over a matrix stored vectorwise by columns
     */
    private class ColMatrixIterator implements Iterator<MatrixEntry> {

        /**
         * Iterates over each column vector
         */
        private SuperIterator<SparseVector, VectorEntry> iterator = new SuperIterator<SparseVector, VectorEntry>(
                Arrays.asList(colD));

        /**
         * Entry returned
         */
        private ColMatrixEntry entry = new ColMatrixEntry();

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
     * Entry of a matrix stored vectorwise by columns
     */
    private static class ColMatrixEntry implements MatrixEntry {

        private int column;

        private VectorEntry entry;

        public void update(int column, VectorEntry entry) {
            this.column = column;
            this.entry = entry;
        }

        public int row() {
            return entry.index();
        }

        public int column() {
            return column;
        }

        public double get() {
            return entry.get();
        }

        public void set(double value) {
            entry.set(value);
        }

    }

}
