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

import java.util.Formatter;
import java.util.Iterator;

/**
 * Partial implementation of <code>Matrix</code>. The following methods throw
 * <code>UnsupportedOperationException</code>, and should be overridden by a
 * subclass:
 * <ul>
 * <li><code>get(int,int)</code></li>
 * <li><code>set(int,int,double)</code></li>
 * <li><code>copy</code></li>
 * <li>All the direct solution methods</li>
 * </ul>
 * <p>
 * For the rest of the methods, simple default implementations using a matrix
 * iterator has been provided. There are some kernel operations which the
 * simpler operations forward to, for instance, <code>mult(Matrix,Matrix)</code>
 * forwards to <code>multAdd(double,Matrix,Matrix)</code>. Subclasses can
 * thus focus on overriding the kernel operations, which are:
 * <ul>
 * <li> <code>multAdd(double,Vector,Vector)</code> and
 * <code>transMultAdd(double,Vector,Vector)</code>. </li>
 * <li> <code>rank1(double,Vector,Vector)</code> and
 * <code>rank1(double,Vector,Vector)</code>.</li>
 * <li> <code>multAdd(double,Matrix,Matrix)</code>,
 * <code>transAmultAdd(double,Matrix,Matrix)</code>,
 * <code>transBmultAdd(double,Matrix,Matrix)</code>, and
 * <code>transABmultAdd(double,Matrix,Matrix)</code>. </li>
 * <li> <code>scale(double)</code>. </li>
 * <li> <code>set(double,Matrix)</code> and <code>add(double,Matrix)</code>.
 * </li>
 * <li> <code>transpose</code> and <code>transpose(Matrix)</code>. </li>
 * <li> All the norms.</li>
 * </ul>
 * <p>
 * Finally, a default iterator is provided by this class, which works by calling
 * the <code>get</code> function. A tailored replacement should be used by
 * subclasses.
 * </ul>
 */
public abstract class AbstractMatrix implements Matrix {

    /**
     * Number of rows
     */
    protected int numRows;

    /**
     * Number of columns
     */
    protected int numColumns;

    /**
     * Constructor for AbstractMatrix
     */
    protected AbstractMatrix(int numRows, int numColumns) {
        if (numRows < 0 || numColumns < 0)
            throw new IndexOutOfBoundsException(
                    "Matrix size cannot be negative");
        this.numRows = numRows;
        this.numColumns = numColumns;
    }

    /**
     * Constructor for AbstractMatrix, same size as A. The invoking constructor
     * should set this matrix equal the argument matrix
     */
    protected AbstractMatrix(Matrix A) {
        this(A.numRows(), A.numColumns());
    }

    public int numRows() {
        return numRows;
    }

    public int numColumns() {
        return numColumns;
    }

    public boolean isSquare() {
        return numRows == numColumns;
    }

    public void set(int row, int column, double value) {
        throw new UnsupportedOperationException();
    }

    public void add(int row, int column, double value) {
        set(row, column, value + get(row, column));
    }

    public double get(int row, int column) {
        throw new UnsupportedOperationException();
    }

    /**
     * Checks the passed row and column indices
     */
    protected void check(int row, int column) {
        if (row < 0)
            throw new IndexOutOfBoundsException("row index is negative (" + row
                    + ")");
        if (column < 0)
            throw new IndexOutOfBoundsException("column index is negative ("
                    + column + ")");
        if (row >= numRows)
            throw new IndexOutOfBoundsException("row index >= numRows (" + row
                    + " >= " + numRows + ")");
        if (column >= numColumns)
            throw new IndexOutOfBoundsException("column index >= numColumns ("
                    + column + " >= " + numColumns + ")");
    }

    public Matrix copy() {
        throw new UnsupportedOperationException();
    }

    public Matrix zero() {
        for (MatrixEntry e : this)
            e.set(0);
        return this;
    }

    public Vector mult(Vector x, Vector y) {
        return mult(1, x, y);
    }

    public Vector mult(double alpha, Vector x, Vector y) {
        return multAdd(alpha, x, y.zero());
    }

    public Vector multAdd(Vector x, Vector y) {
        return multAdd(1, x, y);
    }

    public Vector multAdd(double alpha, Vector x, Vector y) {
        checkMultAdd(x, y);

        if (alpha != 0)
            for (MatrixEntry e : this)
                y.add(e.row(), alpha * e.get() * x.get(e.column()));

        return y;
    }

    /**
     * Checks the arguments to <code>mult</code> and <code>multAdd</code>
     */
    protected void checkMultAdd(Vector x, Vector y) {
        if (numColumns != x.size())
            throw new IndexOutOfBoundsException("A.numColumns != x.size ("
                    + numColumns + " != " + x.size() + ")");
        if (numRows != y.size())
            throw new IndexOutOfBoundsException("A.numRows != y.size ("
                    + numRows + " != " + y.size() + ")");
    }

    public Vector transMult(Vector x, Vector y) {
        return transMult(1, x, y);
    }

    public Vector transMult(double alpha, Vector x, Vector y) {
        return transMultAdd(alpha, x, y.zero());
    }

    public Vector transMultAdd(Vector x, Vector y) {
        return transMultAdd(1, x, y);
    }

    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        checkTransMultAdd(x, y);

        if (alpha != 0)
            for (MatrixEntry e : this)
                y.add(e.column(), alpha * e.get() * x.get(e.row()));

        return y;
    }

    /**
     * Checks the arguments to <code>transMult</code> and
     * <code>transMultAdd</code>
     */
    protected void checkTransMultAdd(Vector x, Vector y) {
        if (numRows != x.size())
            throw new IndexOutOfBoundsException("A.numRows != x.size ("
                    + numRows + " != " + x.size() + ")");
        if (numColumns != y.size())
            throw new IndexOutOfBoundsException("A.numColumns != y.size ("
                    + numColumns + " != " + y.size() + ")");
    }

    public Vector solve(Vector b, Vector x) {
        throw new UnsupportedOperationException();
    }

    public Vector transSolve(Vector b, Vector x) {
        throw new UnsupportedOperationException();
    }

    /**
     * Checks that a matrix inversion is legal for the given arguments. This is
     * for the square case, not for least-squares problems
     */
    protected void checkSolve(Vector b, Vector x) {
        if (!isSquare())
            throw new IndexOutOfBoundsException("!A.isSquare");
        if (numRows != b.size())
            throw new IndexOutOfBoundsException("numRows != b.size (" + numRows
                    + " != " + b.size() + ")");
        if (numColumns != x.size())
            throw new IndexOutOfBoundsException("numColumns != x.size ("
                    + numColumns + " != " + x.size() + ")");
    }

    public Matrix rank1(Vector x) {
        return rank1(1, x);
    }

    public Matrix rank1(double alpha, Vector x) {
        return rank1(alpha, x, x);
    }

    public Matrix rank1(Vector x, Vector y) {
        return rank1(1, x, y);
    }

    public Matrix rank1(double alpha, Vector x, Vector y) {
        checkRank1(x, y);

        if (alpha == 0)
            return this;

        for (VectorEntry ei : x)
            if (ei.get() != 0)
                for (VectorEntry ej : y)
                    if (ej.get() != 0)
                        add(ei.index(), ej.index(), alpha * ei.get() * ej.get());

        return this;
    }

    /**
     * Checks that a vector rank1 update is possible for the given vectors
     */
    protected void checkRank1(Vector x, Vector y) {
        if (!isSquare())
            throw new IndexOutOfBoundsException("!A.isSquare");
        if (x.size() != numRows)
            throw new IndexOutOfBoundsException("x.size != A.numRows ("
                    + x.size() + " != " + numRows + ")");
        if (y.size() != numColumns)
            throw new IndexOutOfBoundsException("y.size != A.numColumns ("
                    + y.size() + " != " + numColumns + ")");
    }

    public Matrix rank2(Vector x, Vector y) {
        return rank2(1, x, y);
    }

    public Matrix rank2(double alpha, Vector x, Vector y) {
        checkRank2(x, y);

        if (alpha == 0)
            return this;

        for (VectorEntry ei : x)
            for (VectorEntry ej : y) {
                add(ei.index(), ej.index(), alpha * ei.get() * ej.get());
                add(ej.index(), ei.index(), alpha * ei.get() * ej.get());
            }

        return this;
    }

    /**
     * Checks that a vector rank2 update is legal with the given vectors
     */
    protected void checkRank2(Vector x, Vector y) {
        if (!isSquare())
            throw new IndexOutOfBoundsException("!A.isSquare");
        if (x.size() != numRows)
            throw new IndexOutOfBoundsException("x.size != A.numRows ("
                    + x.size() + " != " + numRows + ")");
        if (y.size() != numRows)
            throw new IndexOutOfBoundsException("y.size != A.numRows ("
                    + y.size() + " != " + numRows + ")");
    }

    public Matrix mult(Matrix B, Matrix C) {
        return mult(1, B, C);
    }

    public Matrix mult(double alpha, Matrix B, Matrix C) {
        return multAdd(alpha, B, C.zero());
    }

    public Matrix multAdd(Matrix B, Matrix C) {
        return multAdd(1, B, C);
    }

    public Matrix multAdd(double alpha, Matrix B, Matrix C) {
        checkMultAdd(B, C);

        if (alpha != 0)
            for (int i = 0; i < numRows; ++i)
                for (int j = 0; j < C.numColumns(); ++j) {
                    double dot = 0;
                    for (int k = 0; k < numColumns; ++k)
                        dot += get(i, k) * B.get(k, j);
                    C.add(i, j, alpha * dot);
                }

        return C;
    }

    /**
     * Checks the arguments to <code>mult</code> and <code>multAdd</code>
     */
    protected void checkMultAdd(Matrix B, Matrix C) {
        if (numRows != C.numRows())
            throw new IndexOutOfBoundsException("A.numRows != C.numRows ("
                    + numRows + " != " + C.numRows() + ")");
        if (numColumns != B.numRows())
            throw new IndexOutOfBoundsException("A.numColumns != B.numRows ("
                    + numColumns + " != " + B.numRows() + ")");
        if (B.numColumns() != C.numColumns())
            throw new IndexOutOfBoundsException(
                    "B.numColumns != C.numColumns (" + B.numRows() + " != "
                            + C.numColumns() + ")");
    }

    public Matrix transAmult(Matrix B, Matrix C) {
        return transAmult(1, B, C);
    }

    public Matrix transAmult(double alpha, Matrix B, Matrix C) {
        return transAmultAdd(alpha, B, C.zero());
    }

    public Matrix transAmultAdd(Matrix B, Matrix C) {
        return transAmultAdd(1, B, C);
    }

    public Matrix transAmultAdd(double alpha, Matrix B, Matrix C) {
        checkTransAmultAdd(B, C);

        if (alpha != 0)
            for (int i = 0; i < numColumns; ++i)
                for (int j = 0; j < C.numColumns(); ++j) {
                    double dot = 0;
                    for (int k = 0; k < numRows; ++k)
                        dot += get(k, i) * B.get(k, j);
                    C.add(i, j, alpha * dot);
                }

        return C;
    }

    /**
     * Checks the arguments to <code>transAmult</code> and
     * <code>transAmultAdd</code>
     */
    protected void checkTransAmultAdd(Matrix B, Matrix C) {
        if (numRows != B.numRows())
            throw new IndexOutOfBoundsException("A.numRows != B.numRows ("
                    + numRows + " != " + B.numRows() + ")");
        if (numColumns != C.numRows())
            throw new IndexOutOfBoundsException("A.numColumns != C.numRows ("
                    + numColumns + " != " + C.numRows() + ")");
        if (B.numColumns() != C.numColumns())
            throw new IndexOutOfBoundsException(
                    "B.numColumns != C.numColumns (" + B.numColumns() + " != "
                            + C.numColumns() + ")");
    }

    public Matrix transBmult(Matrix B, Matrix C) {
        return transBmult(1, B, C);
    }

    public Matrix transBmult(double alpha, Matrix B, Matrix C) {
        return transBmultAdd(alpha, B, C.zero());
    }

    public Matrix transBmultAdd(Matrix B, Matrix C) {
        return transBmultAdd(1, B, C);
    }

    public Matrix transBmultAdd(double alpha, Matrix B, Matrix C) {
        checkTransBmultAdd(B, C);

        if (alpha != 0)
            for (int i = 0; i < numRows; ++i)
                for (int j = 0; j < C.numColumns(); ++j) {
                    double dot = 0;
                    for (int k = 0; k < numColumns; ++k)
                        dot += get(i, k) * B.get(j, k);
                    C.add(i, j, alpha * dot);
                }

        return C;
    }

    /**
     * Checks the arguments to <code>transBmult</code> and
     * <code>transBmultAdd</code>
     */
    protected void checkTransBmultAdd(Matrix B, Matrix C) {
        if (numColumns != B.numColumns())
            throw new IndexOutOfBoundsException(
                    "A.numColumns != B.numColumns (" + numColumns + " != "
                            + B.numColumns() + ")");
        if (numRows != C.numRows())
            throw new IndexOutOfBoundsException("A.numRows != C.numRows ("
                    + numRows + " != " + C.numRows() + ")");
        if (B.numRows() != C.numColumns())
            throw new IndexOutOfBoundsException("B.numRows != C.numColumns ("
                    + B.numRows() + " != " + C.numColumns() + ")");
    }

    public Matrix transABmult(Matrix B, Matrix C) {
        return transABmult(1, B, C);
    }

    public Matrix transABmult(double alpha, Matrix B, Matrix C) {
        return transABmultAdd(alpha, B, C.zero());
    }

    public Matrix transABmultAdd(Matrix B, Matrix C) {
        return transABmultAdd(1, B, C);
    }

    public Matrix transABmultAdd(double alpha, Matrix B, Matrix C) {
        checkTransABmultAdd(B, C);

        if (alpha != 0)
            for (int i = 0; i < numColumns; ++i)
                for (int j = 0; j < C.numColumns(); ++j) {
                    double dot = 0;
                    for (int k = 0; k < numRows; ++k)
                        dot += get(k, i) * B.get(j, k);
                    C.add(i, j, alpha * dot);
                }

        return C;
    }

    /**
     * Checks the arguments to <code>transABmultAdd</code> and
     * <code>transABmultAdd</code>
     */
    protected void checkTransABmultAdd(Matrix B, Matrix C) {
        if (numRows != B.numColumns())
            throw new IndexOutOfBoundsException("A.numRows != B.numColumns ("
                    + numRows + " != " + B.numColumns() + ")");
        if (numColumns != C.numRows())
            throw new IndexOutOfBoundsException("A.numColumns != C.numRows ("
                    + numColumns + " != " + C.numRows() + ")");
        if (B.numRows() != C.numColumns())
            throw new IndexOutOfBoundsException("B.numRows != C.numColumns ("
                    + B.numRows() + " != " + C.numColumns() + ")");
    }

    public Matrix solve(Matrix B, Matrix X) {    	
        throw new UnsupportedOperationException();
    }

    public Matrix transSolve(Matrix B, Matrix X) {
        throw new UnsupportedOperationException();
    }

    /**
     * Checks that a matrix inversion is legal for the given arguments. This is
     * for the square case, not for least-squares problems
     */
    protected void checkSolve(Matrix B, Matrix X) {
        if (!isSquare())
            throw new IndexOutOfBoundsException("!A.isSquare");
        if (B.numRows() != numRows)
            throw new IndexOutOfBoundsException("B.numRows != A.numRows ("
                    + B.numRows() + " != " + numRows + ")");
        if (B.numColumns() != X.numColumns())
            throw new IndexOutOfBoundsException(
                    "B.numColumns != X.numColumns (" + B.numColumns() + " != "
                            + X.numColumns() + ")");
        if (X.numRows() != numColumns)
            throw new IndexOutOfBoundsException("X.numRows != A.numColumns ("
                    + X.numRows() + " != " + numColumns + ")");
    }

    public Matrix rank1(Matrix C) {
        return rank1(1, C);
    }

    public Matrix rank1(double alpha, Matrix C) {
        checkRank1(C);

        if (alpha == 0)
            return this;

        return C.transBmultAdd(alpha, C, this);
    }

    /**
     * Checks that a matrix rank1 update is possible for the given matrix
     */
    protected void checkRank1(Matrix C) {
        if (!isSquare())
            throw new IndexOutOfBoundsException("!A.isSquare");
        if (numRows != C.numRows())
            throw new IndexOutOfBoundsException("A.numRows != C.numRows ("
                    + numRows + " != " + C.numRows() + ")");
    }

    public Matrix transRank1(Matrix C) {
        return transRank1(1, C);
    }

    public Matrix transRank1(double alpha, Matrix C) {
        checkTransRank1(C);

        if (alpha == 0)
            return this;

        return C.transAmultAdd(alpha, C, this);
    }

    /**
     * Checks that a transposed rank1 update is leagal with the given argument
     */
    protected void checkTransRank1(Matrix C) {
        if (!isSquare())
            throw new IndexOutOfBoundsException("!A.isSquare");
        if (numRows != C.numColumns())
            throw new IndexOutOfBoundsException("A.numRows != C.numColumns ("
                    + numRows + " != " + C.numColumns() + ")");
    }

    public Matrix rank2(Matrix B, Matrix C) {
        return rank2(1, B, C);
    }

    public Matrix rank2(double alpha, Matrix B, Matrix C) {
        checkRank2(B, C);

        if (alpha == 0)
            return this;

        return B.transBmultAdd(alpha, C, C.transBmultAdd(alpha, B, this));
    }

    /**
     * Checks that a rank2 update is legal for the given arguments
     */
    protected void checkRank2(Matrix B, Matrix C) {
        if (!isSquare())
            throw new IndexOutOfBoundsException("!A.isSquare");
        if (B.numRows() != C.numRows())
            throw new IndexOutOfBoundsException("B.numRows != C.numRows ("
                    + B.numRows() + " != " + C.numRows() + ")");
        if (B.numColumns() != C.numColumns())
            throw new IndexOutOfBoundsException(
                    "B.numColumns != C.numColumns (" + B.numColumns() + " != "
                            + C.numColumns() + ")");
    }

    public Matrix transRank2(Matrix B, Matrix C) {
        return transRank2(1, B, C);
    }

    public Matrix transRank2(double alpha, Matrix B, Matrix C) {
        checkTransRank2(B, C);

        if (alpha == 0)
            return this;

        return B.transAmultAdd(alpha, C, C.transAmultAdd(alpha, B, this));
    }

    /**
     * Checks that a transposed rank2 update is leagal with the given arguments
     */
    protected void checkTransRank2(Matrix B, Matrix C) {
        if (!isSquare())
            throw new IndexOutOfBoundsException("!A.isSquare");
        if (numRows != B.numColumns())
            throw new IndexOutOfBoundsException("A.numRows != B.numColumns ("
                    + numRows + " != " + B.numColumns() + ")");
        if (B.numRows() != C.numRows())
            throw new IndexOutOfBoundsException("B.numRows != C.numRows ("
                    + B.numRows() + " != " + C.numRows() + ")");
        if (B.numColumns() != C.numColumns())
            throw new IndexOutOfBoundsException(
                    "B.numColumns != C.numColumns (" + B.numColumns() + " != "
                            + C.numColumns() + ")");
    }

    public Matrix scale(double alpha) {
        if (alpha == 1)
            return this;
        else if (alpha == 0)
            return zero();

        for (MatrixEntry e : this)
            e.set(alpha * e.get());

        return this;
    }

    public Matrix set(Matrix B) {
        return set(1, B);
    }

    public Matrix set(double alpha, Matrix B) {
        checkSize(B);

        if (alpha == 0.)
            return zero();
        if (B == this)
            return scale(alpha);

        zero();
        for (MatrixEntry e : B)
          if (e.get() != 0) set(e.row(), e.column(), alpha * e.get());

        return this;
    }

    public Matrix add(Matrix B) {
        return add(1, B);
    }

    public Matrix add(double alpha, Matrix B) {
        checkSize(B);

        if (alpha != 0)
            for (MatrixEntry e : B)
                add(e.row(), e.column(), alpha * e.get());

        return this;
    }

    /**
     * Checks that the sizes of this matrix and the given conform
     */
    protected void checkSize(Matrix B) {
        if (numRows != B.numRows())
            throw new IndexOutOfBoundsException("A.numRows != B.numRows ("
                    + numRows + " != " + B.numRows() + ")");
        if (numColumns != B.numColumns())
            throw new IndexOutOfBoundsException(
                    "A.numColumns != B.numColumns (" + numColumns + " != "
                            + B.numColumns() + ")");
    }

    public Matrix transpose() {
        checkTranspose();

        for (int j = 0; j < numColumns; ++j)
            for (int i = j + 1; i < numRows; ++i) {
                double value = get(i, j);
                set(i, j, get(j, i));
                set(j, i, value);
            }

        return this;
    }

    /**
     * Checks that the matrix may be transposed
     */
    protected void checkTranspose() {
        if (!isSquare())
            throw new IndexOutOfBoundsException("!A.isSquare");
    }

    public Matrix transpose(Matrix B) {
        checkTranspose(B);

        if (B == this)
            return transpose();

        B.zero();
        for (MatrixEntry e : this)
            B.set(e.column(), e.row(), e.get());

        return B;
    }

    /**
     * Checks that this matrix can be transposed into the given matrix
     */
    protected void checkTranspose(Matrix B) {
        if (numRows != B.numColumns())
            throw new IndexOutOfBoundsException("A.numRows != B.numColumns ("
                    + numRows + " != " + B.numColumns() + ")");
        if (numColumns != B.numRows())
            throw new IndexOutOfBoundsException("A.numColumns != B.numRows ("
                    + numColumns + " != " + B.numRows() + ")");
    }

    public double norm(Norm type) {
        if (type == Norm.One)
            return norm1();
        else if (type == Norm.Frobenius)
            return normF();
        else if (type == Norm.Infinity)
            return normInf();
        else
            // Maxvalue
            return max();
    }

    /**
     * Computes the 1 norm
     */
    protected double norm1() {
        double[] rowSum = new double[numRows];
        for (MatrixEntry e : this)
            rowSum[e.row()] += Math.abs(e.get());
        return max(rowSum);
    }

    /**
     * Computes the Frobenius norm. This implementation is overflow resistant
     */
    protected double normF() {
        double scale = 0, ssq = 1;
        for (MatrixEntry e : this) {
            double Aval = e.get();
            if (Aval != 0) {
                double absxi = Math.abs(Aval);
                if (scale < absxi) {
                    ssq = 1 + ssq * Math.pow(scale / absxi, 2);
                    scale = absxi;
                } else
                    ssq = ssq + Math.pow(absxi / scale, 2);
            }
        }
        return scale * Math.sqrt(ssq);
    }

    /**
     * Computes the infinity norm
     */
    protected double normInf() {
        double[] columnSum = new double[numColumns];
        for (MatrixEntry e : this)
            columnSum[e.column()] += Math.abs(e.get());
        return max(columnSum);
    }

    /**
     * Returns the largest absolute value
     */
    protected double max() {
        double max = 0;
        for (MatrixEntry e : this)
            max = Math.max(Math.abs(e.get()), max);
        return max;
    }

    /**
     * Returns the largest element of the passed array
     */
    protected double max(double[] x) {
        double max = 0;
        for (int i = 0; i < x.length; ++i)
            max = Math.max(x[i], max);
        return max;
    }

    @Override
    public String toString() {
        // Output into coordinate format. Indices start from 1 instead of 0
        Formatter out = new Formatter();

        out.format("%10d %10d %19d\n", numRows, numColumns, Matrices
                .cardinality(this));

        for (MatrixEntry e : this)
            if (e.get() != 0)
                out.format("%10d %10d % .12e\n", e.row() + 1, e.column() + 1, e.get());

        return out.toString();
    }

    public Iterator<MatrixEntry> iterator() {
        return new RefMatrixIterator();
    }

    /**
     * Iterator over a general matrix. Uses column-major traversal
     */
    class RefMatrixIterator implements Iterator<MatrixEntry> {

        /**
         * Matrix cursor
         */
        int row, column;

        /**
         * Matrix entry
         */
        final RefMatrixEntry entry = new RefMatrixEntry();

        public boolean hasNext() {
            return (row < numRows) && (column < numColumns);
        }

        public MatrixEntry next() {
            entry.update(row, column);

            // Traversal first down the columns, then the rows
            if (row < numRows - 1)
                row++;
            else {
                column++;
                row = 0;
            }

            return entry;
        }

        public void remove() {
            entry.set(0);
        }

    }

    /**
     * Matrix entry backed by the matrix. May be reused for higher performance
     */
    class RefMatrixEntry implements MatrixEntry {

        /**
         * Matrix position
         */
        private int row, column;

        /**
         * Updates the entry
         */
        public void update(int row, int column) {
            this.row = row;
            this.column = column;
        }

        public int row() {
            return row;
        }

        public int column() {
            return column;
        }

        public double get() {
            return AbstractMatrix.this.get(row, column);
        }

        public void set(double value) {
            AbstractMatrix.this.set(row, column, value);
        }
    }

}
