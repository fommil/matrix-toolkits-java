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

import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;

import no.uib.cipr.matrix.io.MatrixVectorReader;
import no.uib.cipr.matrix.io.VectorInfo;
import no.uib.cipr.matrix.io.VectorSize;
import no.uib.cipr.matrix.io.VectorInfo.VectorField;

/**
 * Dense vector. Stored by a <code>double[]</code> array of the same length as
 * the vector itself.
 */
public class DenseVector extends AbstractVector implements Serializable {

    /** just the private data */
	private static final long serialVersionUID = 5358813524094629362L;
	
	/**
     * Vector data
     */
    private final double[] data;

    /**
     * Constructor for DenseVector
     * 
     * @param r
     *            Reader to get vector from
     */
    public DenseVector(MatrixVectorReader r) throws IOException {
        // Start with a zero-sized vector
        super(0);

        // Get vector information. Use the header if present, else use a safe
        // default
        VectorInfo info = null;
        if (r.hasInfo())
            info = r.readVectorInfo();
        else
            info = new VectorInfo(true, VectorField.Real);
        VectorSize size = r.readVectorSize(info);

        // Resize the vector to correct size
        this.size = size.size();
        data = new double[size.size()];

        // Check that the vector is in an acceptable format
        if (info.isPattern())
            throw new UnsupportedOperationException(
                    "Pattern vectors are not supported");
        if (info.isComplex())
            throw new UnsupportedOperationException(
                    "Complex vectors are not supported");

        // Read the entries
        if (info.isCoordinate()) {

            // Read coordinate data
            int nz = size.numEntries();
            int[] index = new int[nz];
            double[] entry = new double[nz];
            r.readCoordinate(index, entry);

            // Shift indices from 1-offset to 0-offset
            r.add(-1, index);

            // Store them
            for (int i = 0; i < nz; ++i)
                set(index[i], entry[i]);

        } else
            // info.isArray()
            r.readArray(data);
    }

    /**
     * Constructor for DenseVector
     * 
     * @param size
     *            Size of the vector
     */
    public DenseVector(int size) {
        super(size);
        data = new double[size];
    }

    /**
     * Constructor for DenseVector
     * 
     * @param x
     *            Copies contents from this vector. A deep copy is made
     */
    public DenseVector(Vector x) {
        this(x, true);
    }

    /**
     * Constructor for DenseVector
     * 
     * @param x
     *            Copies contents from this vector
     * @param deep
     *            True for a deep copy. For a shallow copy, <code>x</code>
     *            must be a <code>DenseVector</code>
     */
    public DenseVector(Vector x, boolean deep) {
        super(x);

        if (deep) {
            data = new double[size];
            set(x);
        } else
            data = ((DenseVector) x).getData();
    }

    /**
     * Constructor for DenseVector
     * 
     * @param x
     *            Copies contents from this array
     * @param deep
     *            True for a deep copy. For a shallow copy, <code>x</code> is
     *            aliased with the internal storage
     */
    public DenseVector(double[] x, boolean deep) {
        super(x.length);

        if (deep)
            data = x.clone();
        else
            data = x;
    }

    /**
     * Constructor for DenseVector
     * 
     * @param x
     *            Copies contents from this array in a deep copy
     */
    public DenseVector(double[] x) {
        this(x, true);
    }

    @Override
    public void set(int index, double value) {
        check(index);
        data[index] = value;
    }

    @Override
    public void add(int index, double value) {
        check(index);
        data[index] += value;
    }

    @Override
    public double get(int index) {
        check(index);
        return data[index];
    }

    @Override
    public DenseVector copy() {
        return new DenseVector(this);
    }

    @Override
    public DenseVector zero() {
        Arrays.fill(data, 0);
        return this;
    }

    @Override
    public DenseVector scale(double alpha) {
        for (int i = 0; i < size; ++i)
            data[i] *= alpha;
        return this;
    }

    @Override
    public Vector set(Vector y) {
        if (!(y instanceof DenseVector))
            return super.set(y);

        checkSize(y);

        double[] yd = ((DenseVector) y).getData();
        System.arraycopy(yd, 0, data, 0, size);

        return this;
    }

    @Override
    public Vector set(double alpha, Vector y) {
        if (!(y instanceof DenseVector))
            return super.set(alpha, y);

        checkSize(y);

        if (alpha == 0)
            return zero();

        double[] yd = ((DenseVector) y).getData();

        for (int i = 0; i < size; ++i)
            data[i] = alpha * yd[i];

        return this;
    }

    @Override
    public Vector add(Vector y) {
        if (!(y instanceof DenseVector))
            return super.add(y);

        checkSize(y);

        double[] yd = ((DenseVector) y).getData();

        for (int i = 0; i < size; i++)
            data[i] += yd[i];

        return this;
    }

    @Override
    public Vector add(double alpha, Vector y) {
        if (!(y instanceof DenseVector))
            return super.add(alpha, y);

        checkSize(y);

        if (alpha == 0)
            return this;

        double[] yd = ((DenseVector) y).getData();

        for (int i = 0; i < size; i++)
            data[i] += alpha * yd[i];

        return this;
    }

    @Override
    public double dot(Vector y) {
        if (!(y instanceof DenseVector))
            return super.dot(y);

        checkSize(y);

        double[] yd = ((DenseVector) y).getData();

        double dot = 0.;
        for (int i = 0; i < size; ++i)
            dot += data[i] * yd[i];
        return dot;
    }

    @Override
    protected double norm1() {
        double sum = 0;
        for (int i = 0; i < size; ++i)
            sum += Math.abs(data[i]);
        return sum;
    }

    @Override
    protected double norm2() {
        double norm = 0;
        for (int i = 0; i < size; ++i)
            norm += data[i] * data[i];
        return Math.sqrt(norm);
    }

    @Override
    protected double norm2_robust() {
        double scale = 0, ssq = 1;
        for (int i = 0; i < size; ++i)
            if (data[i] != 0) {
                double absxi = Math.abs(data[i]);
                if (scale < absxi) {
                    ssq = 1 + ssq * (scale / absxi) * (scale / absxi);
                    scale = absxi;
                } else
                    ssq += (absxi / scale) * (absxi / scale);
            }
        return scale * Math.sqrt(ssq);
    }

    @Override
    protected double normInf() {
        double max = 0;
        for (int i = 0; i < size; ++i)
            max = Math.max(Math.abs(data[i]), max);
        return max;
    }

    /**
     * Returns the internal vector contents. The array indices correspond to the
     * vector indices
     */
    public double[] getData() {
        return data;
    }

}
