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
 * Contains the size of a vectir stored in a variant of the <a
 * href="http://math.nist.gov/MatrixMarket">Matrix Market</a> exchange format
 */
public class VectorSize {

    /**
     * Size of the vector
     */
    private int size;

    /**
     * Number of entries stored
     */
    private int numEntries;

    /**
     * Constructor for VectorSize. Assumes dense format
     * 
     * @param size
     *            Size of the matrix
     */
    public VectorSize(int size) {
        this.size = size;
        numEntries = size;

        if (size < 0)
            throw new IllegalArgumentException("size < 0");
    }

    /**
     * Constructor for VectorSize
     * 
     * @param size
     *            Size of the matrix
     * @param numEntries
     *            Number of entries stored
     */
    public VectorSize(int size, int numEntries) {
        this.size = size;
        this.numEntries = numEntries;

        if (size < 0 || numEntries < 0)
            throw new IllegalArgumentException("size < 0 || numEntries < 0");
        if (numEntries > size)
            throw new IllegalArgumentException("numEntries > size");
    }

    /**
     * Returns the size of the vector
     */
    public int size() {
        return size;
    }

    /**
     * Returns the number of entries in the vector
     */
    public int numEntries() {
        return numEntries;
    }
}
