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
 * Partial implementation of a packed matrix
 */
abstract class AbstractPackMatrix extends AbstractMatrix {

    /**
     * Matrix contents
     */
    double[] data;

    /**
     * Matrix is square, so this is either numRows or numColumns
     */
    int n;

    /**
     * Constructor for AbstractPackMatrix
     * 
     * @param n
     *            Size of the matrix. Since the matrix must be square, this
     *            equals both the number of rows and columns
     */
    public AbstractPackMatrix(int n) {
        super(n, n);
        this.n = numRows;
        data = new double[(n * n + n) / 2];
    }

    /**
     * Constructor for AbstractPackMatrix
     * 
     * @param A
     *            Matrix to copy
     */
    public AbstractPackMatrix(Matrix A) {
        this(A, true);
    }

    /**
     * Constructor for AbstractPackMatrix
     * 
     * @param A
     *            A square matrix to copy from
     * @param deep
     *            True for a deep copy, false for a shallow (reference) copy.
     *            References must be instances of a packed matrix
     */
    public AbstractPackMatrix(Matrix A, boolean deep) {
        super(A);

        if (!isSquare())
            throw new IllegalArgumentException("Packed matrix must be square");
        n = A.numRows();
        if (deep) {
            data = new double[(n * n + n) / 2];
            copy(A);
        } else
            this.data = ((AbstractPackMatrix) A).getData();
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
    public Matrix set(Matrix B) {
        if (!(B instanceof AbstractPackMatrix))
            return super.set(B);

        checkSize(B);

        double[] Bd = ((AbstractPackMatrix) B).getData();

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

}
