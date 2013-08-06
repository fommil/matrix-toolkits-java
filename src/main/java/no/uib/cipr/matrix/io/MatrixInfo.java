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
 * Contains information on a matrix in the <a
 * href="http://math.nist.gov/MatrixMarket">Matrix Market</a> exchange format.
 * Supports all valid matrices.
 */
public class MatrixInfo {

    /**
     * What kind of numbers are stored
     */
    public enum MatrixField {

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
     * Symmetry structure of the matrix, if any
     */
    public enum MatrixSymmetry {

        /**
         * General matrix, no symmetry
         */
        General,

        /**
         * Symmetrical matrix
         */
        Symmetric,

        /**
         * Skew symmetrical matrix
         */
        SkewSymmetric,

        /**
         * Hermitian matrix. Only applicable for complex entris
         */
        Hermitian;
    }

    /**
     * True if the matrix is sparse, else false
     */
    private boolean sparse;

    /**
     * Type of data stored
     */
    private MatrixField field;

    /**
     * Matrix symmetry
     */
    private MatrixSymmetry symmetry;

    /**
     * Creates a specific type
     * 
     * @param sparse
     *            True for sparse matrices, else false
     * @param field
     *            Type of data stored
     * @param symmetry
     *            Matrix symmetry
     */
    public MatrixInfo(boolean sparse, MatrixField field, MatrixSymmetry symmetry) {
        this.sparse = sparse;
        this.field = field;
        this.symmetry = symmetry;

        validate();
    }

    /**
     * Validates the representation
     */
    private void validate() {
        if (isDense() && isPattern())
            throw new IllegalArgumentException(
                    "Matrix cannot be dense with pattern storage");
        if (isReal() && isHermitian())
            throw new IllegalArgumentException(
                    "Data cannot be real with hermitian symmetry");
        if (!isComplex() && isHermitian())
            throw new IllegalArgumentException(
                    "Data must be complex with hermitian symmetry");
        if (isPattern() && isSkewSymmetric())
            throw new IllegalArgumentException(
                    "Storage cannot be pattern and skew symmetrical");
    }

    /**
     * Returns <code>true</code> if the matrix is in coordinate format, else
     * <code>false</code>
     */
    public boolean isSparse() {
        return sparse;
    }

    /**
     * Returns <code>true</code> if the matrix is in coordinate format, else
     * <code>false</code>
     */
    public boolean isCoordinate() {
        return sparse;
    }

    /**
     * Returns <code>true</code> if the matrix is in array format, else
     * <code>false</code>
     */
    public boolean isDense() {
        return !sparse;
    }

    /**
     * Returns <code>true</code> if the matrix is in array format, else
     * <code>false</code>
     */
    public boolean isArray() {
        return !sparse;
    }

    /**
     * Returns <code>true</code> if the matrix stores real numbers, else
     * <code>false</code>
     */
    public boolean isReal() {
        return field == MatrixField.Real;
    }

    /**
     * Returns <code>true</code> if the matrix stores integers, else
     * <code>false</code>
     */
    public boolean isInteger() {
        return field == MatrixField.Integer;
    }

    /**
     * Returns <code>true</code> if the matrix stores complex numbers, else
     * <code>false</code>
     */
    public boolean isComplex() {
        return field == MatrixField.Complex;
    }

    /**
     * Returns <code>true</code> if the matrix does not store any numbers,
     * else <code>false</code>
     */
    public boolean isPattern() {
        return field == MatrixField.Pattern;
    }

    /**
     * Returns <code>true</code> if the matrix form is general, else
     * <code>false</code>
     */
    public boolean isGeneral() {
        return symmetry == MatrixSymmetry.General;
    }

    /**
     * Returns <code>true</code> if the matrix is symmetrical, else
     * <code>false</code>
     */
    public boolean isSymmetric() {
        return symmetry == MatrixSymmetry.Symmetric;
    }

    /**
     * Returns <code>true</code> if the matrix is skew-symmetrical, else
     * <code>false</code>
     */
    public boolean isSkewSymmetric() {
        return symmetry == MatrixSymmetry.SkewSymmetric;
    }

    /**
     * Returns <code>true</code> if the matrix is Hermitian, else
     * <code>false</code>
     */
    public boolean isHermitian() {
        return symmetry == MatrixSymmetry.Hermitian;
    }

    /**
     * Returns a string representation of the specifier. Can be used to provide
     * a header for writing to a file. It is a two-line output, which can look
     * like this:
     * 
     * <pre>
     *       %%MatrixMarket matrix coordinate real general
     * </pre>
     */
    @Override
    public String toString() {
        StringBuilder buf = new StringBuilder();

        buf.append("%%MatrixMarket matrix ");

        if (isSparse())
            buf.append("coordinate ");
        else
            buf.append("array ");

        if (isReal())
            buf.append("real ");
        else if (isComplex())
            buf.append("complex ");
        else if (isPattern())
            buf.append("pattern ");
        else if (isInteger())
            buf.append("integer ");
        else
            // This should never happen
            throw new IllegalArgumentException("Unknown field specification");

        if (isGeneral())
            buf.append("general\n");
        else if (isSymmetric())
            buf.append("symmetric\n");
        else if (isSkewSymmetric())
            buf.append("skew-symmetric\n");
        else if (isHermitian())
            buf.append("Hermitian\n");
        else
            // This should never happen
            throw new IllegalArgumentException("Unknown symmetry specification");

        return buf.toString();
    }

}
