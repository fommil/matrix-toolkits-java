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

/**
 * Static utility methods for matrices and vectors
 */
public final class Matrices {

    private Matrices() {
        // No need to instantiate
    }

	/**
	 * <code>max(1, M)</code> provided as a convenience for 'leading dimension' calculations.
	 * 
	 * @param n
	 */
	static int ld(int n) {
		return Math.max(1, n);
	}

	/**
	 * <code>max(1, max(M, N))</code> provided as a convenience for 'leading dimension'
	 * calculations.
	 * 
	 * @param m
	 * @param n
	 */
	static int ld(int m, int n) {
		return Math.max(1, Math.max(m, n));
	}
    
    /**
     * Returns the number of non-zero entries in the given vector
     */
    public static int cardinality(Vector x) {
        int nz = 0;
        for (VectorEntry e : x)
            if (e.get() != 0)
                nz++;
        return nz;
    }

    /**
     * Returns the number of non-zero entries in the given matrix
     */
    public static int cardinality(Matrix A) {
        int nz = 0;
        for (MatrixEntry e : A)
            if (e.get() != 0)
                nz++;
        return nz;
    }

    /**
     * Returns an array of arrays containing a copy of the given matrix. Each
     * array contains one row.
     */
    public static double[][] getArray(Matrix A) {
        double[][] Ad = new double[A.numRows()][A.numColumns()];
        for (MatrixEntry e : A)
            Ad[e.row()][e.column()] = e.get();
        return Ad;
    }

    /**
     * Returns a dense array containing a copy of the given vector
     */
    public static double[] getArray(Vector x) {
        double[] xd = new double[x.size()];
        for (VectorEntry e : x)
            xd[e.index()] = e.get();
        return xd;
    }

    /**
     * Returns the identity matrix of the given size
     * 
     * @param size
     *            Number of rows/columns of the matrix
     * @return Matrix of the given size, with ones on the main diagonal
     */
    public static DenseMatrix identity(int size) {
        DenseMatrix A = new DenseMatrix(size, size);
        for (int i = 0; i < size; ++i)
            A.set(i, i, 1);
        return A;
    }

    /**
     * Creates a random vector. Numbers are drawn from a uniform distribution
     * between 0 and 1
     * 
     * @param size
     *            Size of the vector
     */
    public static Vector random(int size) {
        return random(new DenseVector(size));
    }

    /**
     * Populates a vector with random numbers drawn from a uniform distribution
     * between 0 and 1
     * 
     * @param x
     *            Vector to populate
     */
    public static Vector random(Vector x) {
        for (int i = 0; i < x.size(); ++i)
            x.set(i, Math.random());
        return x;
    }

    /**
     * Creates a random matrix. Numbers are drawn from a uniform distribution
     * between 0 and 1
     * 
     * @param numRows
     *            Number of rows
     * @param numColumns
     *            Number of columns
     */
    public static Matrix random(int numRows, int numColumns) {
        return random(new DenseMatrix(numRows, numColumns));
    }

    /**
     * Populates a matrix with random numbers drawn from a uniform distribution
     * between 0 and 1
     * 
     * @param A
     *            Matrix to populate
     */
    public static Matrix random(Matrix A) {
        for (int j = 0; j < A.numColumns(); ++j)
            for (int i = 0; i < A.numRows(); ++i)
                A.set(i, j, Math.random());
        return A;
    }

    /**
     * Returns a synchronized vector which wraps the given vector. Only the
     * <code>set(int, double)</code> and <code>add(int, double)</code>
     * methods and their blocked versions are synchronized.
     * <p>
     * <b>Note: </b> Do not use the wrapped vector for any operations besides
     * matrix assembly, as these operations may be very slow.
     * 
     * @param x
     *            Vector to be wrapped
     * @return A thin wrapper around <code>x</code>
     */
    public static Vector synchronizedVector(Vector x) {
        return new SynchronizedVector(x);
    }

    /**
     * Returns a synchronized matrix which wraps the given matrix. Only the
     * <code>set(int, int, double)</code> and
     * <code>add(int, int, double)</code> methods and their blocked versions
     * are synchronized.
     * <p>
     * <b>Note: </b> Do not use the wrapped matrix for any operations besides
     * matrix assembly, as these operations may be very slow.
     * 
     * @param A
     *            Matrix to be wrapped
     * @return A thin wrapper around <code>A</code>
     */
    public static Matrix synchronizedMatrix(Matrix A) {
        return new SynchronizedMatrix(A);
    }

    /**
     * Returns a synchronized matrix which wraps the given matrix. Only the
     * <code>set(int, int, double)</code> and
     * <code>add(int, int, double)</code> methods and their blocked versions
     * are synchronized.
     * <p>
     * The locking provided is finer than the locking of the whole matrix, as
     * different threads can access different rows simultaneous, while only one
     * thread can access a given row at a time. Use this for row-major matrices,
     * <i>not </i> for column-major matrices.
     * <p>
     * <b>Note: </b> Do not use the wrapped matrix for any operations besides
     * matrix assembly, as these operations may be very slow.
     * 
     * @param A
     *            Matrix to be wrapped
     * @return A thin wrapper around <code>A</code>. Individual rows are
     *         locked
     */
    public static Matrix synchronizedMatrixByRows(Matrix A) {
        return new SynchronizedRowMatrix(A);
    }

    /**
     * Returns a synchronized matrix which wraps the given matrix. Only the
     * <code>set(int, int, double)</code> and
     * <code>add(int, int, double)</code> methods and their blocked versions
     * are synchronized.
     * <p>
     * The locking provided is finer than the locking of the whole matrix, as
     * different threads can access different columns simultaneous, while only
     * one thread can access a given column at a time. Use this for column-major
     * matrices, <i>not </i> for row-major matrices.
     * <p>
     * <b>Note: </b> Do not use the wrapped matrix for any operations besides
     * matrix assembly, as these operations may be very slow.
     * 
     * @param A
     *            Matrix to be wrapped
     * @return A thin wrapper around <code>A</code>. Individual columns are
     *         locked
     */
    public static Matrix synchronizedMatrixByColumns(Matrix A) {
        return new SynchronizedColumnMatrix(A);
    }

    /**
     * Returns a view into the given matrix. This view is only for easing some
     * matrix-assembly cases, not for general use. To extract a more
     * higher-performing and general matrix, create a copy of the submatrix. The
     * result is a {@link no.uib.cipr.matrix.DenseMatrix DenseMatrix}.
     * 
     * @param A
     *            Matrix to create view on
     * @param row
     *            Rows to access. Must be within the bounds of <code>A</code>
     * @param column
     *            Columns to access. Must be within the bounds of <code>A</code>
     * @return Submatrix of <code>A</code>. Changing it will change the
     *         backing matrix
     */
    public static Matrix getSubMatrix(Matrix A, int[] row, int[] column) {
        return new RefMatrix(A, row, column);
    }

    /**
     * Returns a view into the given vector. This view is only for easing some
     * vector-assembly cases, not for general use. To extract a more
     * higher-performing and general vector, create a copy of the subvector. The
     * result is a {@link no.uib.cipr.matrix.DenseVector DenseVector}.
     * 
     * @param x
     *            Vector to create view on
     * @param index
     *            Indices to access. Must be within the bounds of <code>x</code>
     * @return Submatrix of <code>x</code>. Changing it will change the
     *         backing matrix
     */
    public static Vector getSubVector(Vector x, int[] index) {
        return new RefVector(x, index);
    }

    /**
     * Matrix backed by another matrix. Used by <code>getSubMatrix</code>
     */
    private static class RefMatrix extends AbstractMatrix {

        private Matrix A;

        private int[] row, column;

        public RefMatrix(Matrix A, int[] row, int[] column) {
            super(row.length, column.length);
            this.A = A;
            this.row = row;
            this.column = column;
        }

        @Override
        public void add(int row, int column, double value) {
            A.add(this.row[row], this.column[column], value);
        }

        @Override
        public DenseMatrix copy() {
            return new DenseMatrix(this);
        }

        @Override
        public double get(int row, int column) {
            return A.get(this.row[row], this.column[column]);
        }

        @Override
        public void set(int row, int column, double value) {
            A.set(this.row[row], this.column[column], value);
        }

    }

    /**
     * Vector backed by another vector. Used by <code>getSubVector</code>
     */
    private static class RefVector extends AbstractVector {

        private Vector x;

        private int[] index;

        public RefVector(Vector x, int[] index) {
            super(index.length);
            this.x = x;
            this.index = index;
        }

        @Override
        public void add(int index, double value) {
            x.add(this.index[index], value);
        }

        @Override
        public DenseVector copy() {
            return new DenseVector(this);
        }

        @Override
        public double get(int index) {
            return x.get(this.index[index]);
        }

        @Override
        public void set(int index, double value) {
            x.set(this.index[index], value);
        }

    }

    /**
     * Ensures correctness in the vector assembly. Since it extends the
     * AbstractVector class, algebraic operations will be slow. It is not
     * possible to implement Vector and delegate calls to the imbedded vector,
     * since casting to the imbedded vector is not possible
     */
    private static class SynchronizedVector extends AbstractVector {

        private Vector x;

        public SynchronizedVector(Vector x) {
            super(x);
            this.x = x;
        }

        @Override
        public synchronized void add(int index, double value) {
            x.add(index, value);
        }

        @Override
        public synchronized void set(int index, double value) {
            x.set(index, value);
        }

        @Override
        public synchronized double get(int index) {
            return x.get(index);
        }

        @Override
        public Vector copy() {
            return Matrices.synchronizedVector(x.copy());
        }

    }

    /**
     * Ensures correctness in the matrix assembly. Since it extends the
     * AbstractMatrix class, algebraic operations will be slow. It is not
     * possible to implement Matrix and delegate calls to the imbedded matrix,
     * since casting to the imbedded matrix is not possible
     */
    private static class SynchronizedMatrix extends AbstractMatrix {

        private Matrix A;

        public SynchronizedMatrix(Matrix A) {
            super(A);
            this.A = A;
        }

        @Override
        public synchronized void add(int row, int column, double value) {
            A.add(row, column, value);
        }

        @Override
        public synchronized void set(int row, int column, double value) {
            A.set(row, column, value);
        }

        @Override
        public synchronized double get(int row, int column) {
            return A.get(row, column);
        }

        @Override
        public Matrix copy() {
            return Matrices.synchronizedMatrix(A.copy());
        }

    }

    /**
     * Ensures correctness in the matrix assembly. Since it extends the
     * AbstractMatrix class, algebraic operations will be slow. It is not
     * possible to implement Matrix and delegate calls to the imbedded matrix,
     * since casting to the imbedded matrix is not possible
     * <p>
     * Locks individual rows instead of the whole matrix
     */
    private static class SynchronizedRowMatrix extends AbstractMatrix {

        private Matrix A;

        private Object[] lock;

        public SynchronizedRowMatrix(Matrix A) {
            super(A);
            this.A = A;
            lock = new Object[A.numRows()];
            for (int i = 0; i < lock.length; ++i)
                lock[i] = new Object();
        }

        @Override
        public void add(int row, int column, double value) {
            synchronized (lock[row]) {
                A.add(row, column, value);
            }
        }

        @Override
        public void set(int row, int column, double value) {
            synchronized (lock[row]) {
                A.set(row, column, value);
            }
        }

        @Override
        public double get(int row, int column) {
            return A.get(row, column);
        }

        @Override
        public Matrix copy() {
            return Matrices.synchronizedMatrixByRows(A.copy());
        }

    }

    /**
     * Ensures correctness in the matrix assembly. Implements matrix instead of
     * subclassing the abstract matrix in order to correctly delegate every
     * method to possbly overridden method in the encapsulated matrix.
     * <p>
     * Locks individual columns instead of the whole matrix
     */
    private static class SynchronizedColumnMatrix extends AbstractMatrix {

        private Matrix A;

        private Object[] lock;

        public SynchronizedColumnMatrix(Matrix A) {
            super(A);
            this.A = A;
            lock = new Object[A.numColumns()];
            for (int i = 0; i < lock.length; ++i)
                lock[i] = new Object();
        }

        @Override
        public void add(int row, int column, double value) {
            synchronized (lock[column]) {
                A.add(row, column, value);
            }
        }

        @Override
        public void set(int row, int column, double value) {
            synchronized (lock[column]) {
                A.set(row, column, value);
            }
        }

        @Override
        public double get(int row, int column) {
            return A.get(row, column);
        }

        @Override
        public Matrix copy() {
            return Matrices.synchronizedMatrixByColumns(A.copy());
        }

    }

    /**
     * Creates a continuous linear index.
     * 
     * @param from
     *            Start, inclusive
     * @param to
     *            Stop, exclusive
     */
    public static int[] index(int from, int to) {
        int length = to - from;

        if (length < 0)
            length = 0;

        int[] index = new int[length];
        for (int i = from, j = 0; j < length; ++i, ++j)
            index[j] = i;
        return index;
    }

    /**
     * Creates a strided linear index.
     * 
     * @param from
     *            Start, inclusive
     * @param stride
     *            <code>stride=1</code> for continuous. Negative strides are
     *            allowed
     * @param to
     *            Stop, exclusive
     */
    public static int[] index(int from, int stride, int to) {
        if (stride == 1)
            return index(from, to);
        else if (stride == 0)
            return new int[0];

        if (to <= from && stride > 0)
            return new int[0];
        if (from <= to && stride < 0)
            return new int[0];

        int length = Math.abs((to - from) / stride);
        if (Math.abs((to - from) % stride) > 0)
            length++;

        if (length < 0)
            length = 0;

        int[] index = new int[length];
        for (int i = from, j = 0; j < length; i += stride, ++j)
            index[j] = i;
        return index;
    }

    /**
     * Finds the number of non-zero entries on each row
     */
    public static int[] rowBandwidth(Matrix A) {
        int[] nz = new int[A.numRows()];

        for (MatrixEntry e : A)
            nz[e.row()]++;

        return nz;
    }

    /**
     * Finds the number of non-zero entries on each column
     */
    public static int[] columnBandwidth(Matrix A) {
        int[] nz = new int[A.numColumns()];

        for (MatrixEntry e : A)
            nz[e.column()]++;

        return nz;
    }

    /**
     * Finds the number of diagonals below the main diagonal. Useful for
     * converting a general matrix into a banded matrix
     */
    public static int getNumSubDiagonals(Matrix A) {
        int kl = 0;

        for (MatrixEntry e : A)
            kl = Math.max(kl, e.row() - e.column());

        return kl;
    }

    /**
     * Finds the number of diagonals above the main diagonal. Useful for
     * converting a general matrix into a banded matrix
     */
    public static int getNumSuperDiagonals(Matrix A) {
        int ku = 0;

        for (MatrixEntry e : A)
            ku = Math.max(ku, e.column() - e.row());

        return ku;
    }

    /**
     * Sets the selected rows of <code>A</code> equal zero, and puts
     * <code>diagonal</code> on the diagonal of those rows. Useful for
     * enforcing boundary conditions
     */
    public static void zeroRows(Matrix A, double diagonal, int... row) {
        // Sort the rows
        int[] rowS = row.clone();
        Arrays.sort(rowS);

        for (MatrixEntry e : A) {
            int j = java.util.Arrays.binarySearch(rowS, e.row());
            if (j >= 0) { // Found
                if (e.row() == e.column()) // Diagonal
                    e.set(diagonal);
                else
                    // Off diagonal
                    e.set(0);
            }
        }

        // Ensure the diagonal is set. This is necessary in case of missing
        // rows
        if (diagonal != 0)
            for (int rowI : row)
                A.set(rowI, rowI, diagonal);
    }

    /**
     * Sets the selected columns of <code>A</code> equal zero, and puts
     * <code>diagonal</code> on the diagonal of those columns. Useful for
     * enforcing boundary conditions
     */
    public static void zeroColumns(Matrix A, double diagonal, int... column) {
        // Sort the columns
        int[] columnS = column.clone();
        Arrays.sort(columnS);

        for (MatrixEntry e : A) {
            int j = java.util.Arrays.binarySearch(columnS, e.column());
            if (j >= 0) { // Found
                if (e.row() == e.column()) // Diagonal
                    e.set(diagonal);
                else
                    // Off diagonal
                    e.set(0);
            }
        }

        // Ensure the diagonal is set. This is necessary in case of missing
        // columns
        if (diagonal != 0)
            for (int columnI : column)
                A.set(columnI, columnI, diagonal);
    }

  public static DenseVector getColumn(Matrix m, int j) {
    DenseVector v = new DenseVector(m.numRows());
    for (int i = 0; i < v.size(); i++) {
      v.set(i, m.get(i, j));
    }
    return v;
  }

}
