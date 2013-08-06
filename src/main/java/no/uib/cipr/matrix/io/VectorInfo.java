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
 * Contains information on a vector in a variant of the <a
 * href="http://math.nist.gov/MatrixMarket">Matrix Market</a> exchange format
 */
public class VectorInfo {

    /**
     * What kind of numbers are stored
     */
    public enum VectorField {
        /**
         * Real numbers
         */
        Real,

        /**
         * Integers
         */
        Integer,

        /**
         * Complex numbers
         */
        Complex,

        /**
         * Pattern matrix. No numbers stored
         */
        Pattern;
    }

    /**
     * True if the vector is sparse, else false
     */
    private boolean sparse;

    /**
     * Type of data stored
     */
    private VectorField field;

    /**
     * Creates a specific type
     * 
     * @param sparse
     *            True for sparse vectors, else false
     * @param field
     *            Type of data stored
     */
    public VectorInfo(boolean sparse, VectorField field) {
        this.sparse = sparse;
        this.field = field;

        validate();
    }

    /**
     * Validates the representation
     */
    private void validate() {
        if (isDense() && isPattern())
            throw new IllegalArgumentException(
                    "Vector cannot be dense with pattern storage");
    }

    /**
     * Returns <code>true</code> if the vector is in coordinate format, else
     * <code>false</code>
     */
    public boolean isSparse() {
        return sparse;
    }

    /**
     * Returns <code>true</code> if the vector is in coordinate format, else
     * <code>false</code>
     */
    public boolean isCoordinate() {
        return sparse;
    }

    /**
     * Returns <code>true</code> if the vector is in array format, else
     * <code>false</code>
     */
    public boolean isDense() {
        return !sparse;
    }

    /**
     * Returns <code>true</code> if the vector is in array format, else
     * <code>false</code>
     */
    public boolean isArray() {
        return !sparse;
    }

    /**
     * Returns <code>true</code> if the vector stores real numbers, else
     * <code>false</code>
     */
    public boolean isReal() {
        return field == VectorField.Real;
    }

    /**
     * Returns <code>true</code> if the vector stores integers, else
     * <code>false</code>
     */
    public boolean isInteger() {
        return field == VectorField.Integer;
    }

    /**
     * Returns <code>true</code> if the vector stores complex numbers, else
     * <code>false</code>
     */
    public boolean isComplex() {
        return field == VectorField.Complex;
    }

    /**
     * Returns <code>true</code> if the vector does not store any numbers,
     * else <code>false</code>
     */
    public boolean isPattern() {
        return field == VectorField.Pattern;
    }

    /**
     * Returns a string representation of the specifier. Can be used to provide
     * a header for writing to a file. It is a two-line output, which can look
     * like this:
     * 
     * <pre>
     *      %%MatrixMarket vector coordinate real
     * </pre>
     */
    @Override
    public String toString() {
        StringBuilder buf = new StringBuilder();

        buf.append("%%MatrixMarket vector ");

        if (isSparse())
            buf.append("coordinate ");
        else
            buf.append("array ");

        if (isReal())
            buf.append("real\n");
        else if (isComplex())
            buf.append("complex\n");
        else if (isPattern())
            buf.append("pattern\n");
        else if (isInteger())
            buf.append("integer\n");
        else
            // This should never happen
            throw new IllegalArgumentException("Unknown field specification");

        return buf.toString();
    }

}
