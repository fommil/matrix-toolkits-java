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

package no.uib.cipr.matrix.io;

/**
 * Contains the size of a matrix stored in the <a
 * href="http://math.nist.gov/MatrixMarket">Matrix Market</a> exchange format
 */
public class MatrixSize {

    /**
     * Number of rows
     */
    private int numRows;

    /**
     * Number of columns
     */
    private int numColumns;

    /**
     * Number of entries stored
     */
    private int numEntries;

    /**
     * Constructor for MatrixSize
     * 
     * @param numRows
     *            Number of rows in the matrix
     * @param numColumns
     *            Number of columns in the matrix
     * @param info
     *            Info on the matrix
     */
    public MatrixSize(int numRows, int numColumns, MatrixInfo info) {
        this.numRows = numRows;
        this.numColumns = numColumns;

        if (!info.isDense())
            throw new IllegalArgumentException("Matrix must be dense");

        if (info.isGeneral())
            numEntries = numRows * numColumns;
        else if (info.isSymmetric() || info.isHermitian())
            numEntries = ((numRows * numColumns - numRows) / 2 + numRows);
        else if (info.isSkewSymmetric())
            numEntries = (numRows * numColumns - numRows) / 2;
    }

    /**
     * Constructor for MatrixSize
     * 
     * @param numRows
     *            Number of rows in the matrix
     * @param numColumns
     *            Number of columns in the matrix
     * @param numEntries
     *            Number of entries stored
     */
    public MatrixSize(int numRows, int numColumns, int numEntries) {
        this.numRows = numRows;
        this.numColumns = numColumns;
        this.numEntries = numEntries;

        // We do this to avoid overflows
        long maxR = numRows, maxC = numColumns, max = maxR * maxC;
        if (numEntries > max)
            throw new IllegalArgumentException(
                    "numEntries > numRows*numColumns");
    }

    /**
     * Returns the number of rows in the matrix
     */
    public int numRows() {
        return numRows;
    }

    /**
     * Returns the number of columns in the matrix
     */
    public int numColumns() {
        return numColumns;
    }

    /**
     * Returns the number of entries stored
     */
    public int numEntries() {
        return numEntries;
    }

    /**
     * Returns <code>true</code> if the matrix is square, else
     * <code>false</code>
     */
    public boolean isSquare() {
        return numRows == numColumns;
    }

}
