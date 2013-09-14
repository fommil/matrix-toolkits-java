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

import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * Partial implementation of a dense matrix
 */
abstract class AbstractDenseMatrix extends AbstractMatrix {

    /**
     * Matrix contents
     */
    double[] data;

    /**
     * Constructor for AbstractDenseMatrix. The matrix contents will be set to
     * zero
     * 
     * @param numRows
     *            Number of rows
     * @param numColumns
     *            Number of columns
     */
    public AbstractDenseMatrix(int numRows, int numColumns) {
        super(numRows, numColumns);

        // We know that numRows and numColumns is positive from the super
        // constructor.
        
        final long size = (long) numRows * numColumns;
        if (size > Integer.MAX_VALUE)
        {
            throw new IllegalArgumentException("Matrix of "
                + numRows + " x " + numColumns
                + " = " + size + " elements is too large to be allocated using a single Java array.");
        }

        data = new double[numRows * numColumns];
    }

    /**
     * Constructor for AbstractDenseMatrix. Matrix is copied from the supplied
     * matrix
     * 
     * @param A
     *            Matrix to copy from
     */
    public AbstractDenseMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for AbstractDenseMatrix. Matrix is copied from the supplied
     * matrix
     * 
     * @param A
     *            Matrix to copy from
     * @param deep
     *            True for deep copy, false for reference
     */
    public AbstractDenseMatrix(Matrix A, boolean deep) {
        super(A);

        if (deep) {
            data = new double[numRows * numColumns];
            copy(A);
        } else
            this.data = ((AbstractDenseMatrix) A).getData();
    }

    /**
     * Set this matrix equal to the given matrix
     */
    abstract void copy(Matrix A);

    /**
     * Returns the matrix contents. Ordering depends on the underlying storage
     * assumptions
     */
    public double[] getData() {
        return data;
    }

    @Override
    public void add(int row, int column, double value) {
        data[getIndex(row, column)] += value;
    }

    @Override
    public void set(int row, int column, double value) {
        data[getIndex(row, column)] = value;
    }

    @Override
    public double get(int row, int column) {
        return data[getIndex(row, column)];
    }

    /**
     * Checks the row and column indices, and returns the linear data index
     */
    int getIndex(int row, int column) {
        check(row, column);
        return row + column * numRows;
    }

    @Override
    public Matrix set(Matrix B) {
        // using instanceof results in weird problems
        // with implementations that mask some values
        if (!(getClass().isAssignableFrom(B.getClass())))
            return super.set(B);

        checkSize(B);

        double[] Bd = ((AbstractDenseMatrix) B).getData();

        if (Bd == data)
            return this;

        System.arraycopy(Bd, 0, data, 0, data.length);

        return this;
    }

    @Override
    public Matrix zero() {
        Arrays.fill(data, 0);
        return this;
    }

  @Override
  public String toString() {
    StringBuilder out = new StringBuilder();
    DecimalFormat df = new DecimalFormat("####0.00");

    for (int i = 0; i < numRows(); i++) {
      for (int j = 0; j < numColumns(); j++) {
        double value = get(i, j);
        if (value >= 0) out.append(" ");
        out.append(" " + df.format(value));
      }
      out.append("\n");
    }

    return out.toString();
  }

}
