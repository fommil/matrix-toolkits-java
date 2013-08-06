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

import java.io.Serializable;
import java.util.Formatter;
import java.util.Iterator;

/**
 * Partial implementation of <code>Vector</code>. The following methods throw
 * <code>UnsupportedOperationException</code>, and should be overridden by a
 * subclass:
 * <ul>
 * <li><code>get(int)</code></li>
 * <li><code>set(int,double)</code></li>
 * <li><code>copy</code></li>
 * </ul>
 * <p>
 * For the rest of the methods, simple default implementations using a vector
 * iterator has been provided. There are some kernel operations which the
 * simpler operations forward to, and they are:
 * <ul>
 * <li> <code>add(double,Vector)</code> and <code>set(double,Vector)</code>.
 * </li>
 * <li> <code>scale(double)</code>.</li>
 * <li><code>dot(Vector)</code> and all the norms. </li>
 * </ul>
 * <p>
 * Finally, a default iterator is provided by this class, which works by calling
 * the <code>get</code> function. A tailored replacement should be used by
 * subclasses.
 * </ul>
 */
public abstract class AbstractVector implements Vector, Serializable {

    /**
     * Size of the vector
     */
    protected int size;

    /**
     * Constructor for AbstractVector.
     * 
     * @param size
     *            Size of the vector
     */
    protected AbstractVector(int size) {
        if (size < 0)
            throw new IllegalArgumentException("Vector size cannot be negative");
        this.size = size;
    }

    /**
     * Constructor for AbstractVector, same size as x
     * 
     * @param x
     *            Vector to get the size from
     */
    protected AbstractVector(Vector x) {
        this.size = x.size();
    }

    public int size() {
        return size;
    }

    public void set(int index, double value) {
        throw new UnsupportedOperationException();
    }

    public void add(int index, double value) {
        set(index, value + get(index));
    }

    public double get(int index) {
        throw new UnsupportedOperationException();
    }

    public Vector copy() {
        throw new UnsupportedOperationException();
    }

    /**
     * Checks the index
     */
    protected void check(int index) {
        if (index < 0)
            throw new IndexOutOfBoundsException("index is negative (" + index
                    + ")");
        if (index >= size)
            throw new IndexOutOfBoundsException("index >= size (" + index
                    + " >= " + size + ")");
    }

    public Vector zero() {
        for (VectorEntry e : this)
            e.set(0);
        return this;
    }

    public Vector scale(double alpha) {
        if (alpha == 0)
            return zero();
        else if (alpha == 1)
            return this;

        for (VectorEntry e : this)
            e.set(alpha * e.get());

        return this;
    }

    public Vector set(Vector y) {
        return set(1, y);
    }

    public Vector set(double alpha, Vector y) {
        checkSize(y);

        if (alpha == 0)
            return zero();

        zero();
        for (VectorEntry e : y)
            set(e.index(), alpha * e.get());

        return this;
    }

    public Vector add(Vector y) {
        return add(1, y);
    }

    public Vector add(double alpha, Vector y) {
        checkSize(y);

        if (alpha == 0)
            return this;

        for (VectorEntry e : y)
            add(e.index(), alpha * e.get());

        return this;
    }

    public double dot(Vector y) {
        checkSize(y);

        double ret = 0;
        for (VectorEntry e : this)
            ret += e.get() * y.get(e.index());
        return ret;
    }

    /**
     * Checks for conformant sizes
     */
    protected void checkSize(Vector y) {
        if (size != y.size())
            throw new IndexOutOfBoundsException("x.size != y.size (" + size
                    + " != " + y.size() + ")");
    }

    public double norm(Norm type) {
        if (type == Norm.One)
            return norm1();
        else if (type == Norm.Two)
            return norm2();
        else if (type == Norm.TwoRobust)
            return norm2_robust();
        else
            // Infinity
            return normInf();
    }

    protected double norm1() {
        double sum = 0;
        for (VectorEntry e : this)
            sum += Math.abs(e.get());
        return sum;
    }

    protected double norm2() {
        double norm = 0;
        for (VectorEntry e : this)
            norm += e.get() * e.get();
        return Math.sqrt(norm);
    }

    protected double norm2_robust() {
        double scale = 0, ssq = 1;
        for (VectorEntry e : this) {
            double xval = e.get();
            if (xval != 0) {
                double absxi = Math.abs(xval);
                if (scale < absxi) {
                    ssq = 1 + ssq * Math.pow(scale / absxi, 2);
                    scale = absxi;
                } else
                    ssq = ssq + Math.pow(absxi / scale, 2);
            }
        }
        return scale * Math.sqrt(ssq);
    }

    protected double normInf() {
        double max = 0;
        for (VectorEntry e : this)
            max = Math.max(Math.abs(e.get()), max);
        return max;
    }

    public Iterator<VectorEntry> iterator() {
        return new RefVectorIterator();
    }

    @Override
    public String toString() {
        // Output into coordinate format. Indices start from 1 instead of 0
        Formatter out = new Formatter();

        out.format("%10d %19d\n", size, Matrices.cardinality(this));

        for (VectorEntry e : this)
            if (e.get() != 0)
                out.format("%10d % .12e\n", e.index() + 1, e.get());

        return out.toString();
    }

    /**
     * Iterator over a general vector
     */
    private class RefVectorIterator implements Iterator<VectorEntry> {

        private int index;

        private final RefVectorEntry entry = new RefVectorEntry();

        public boolean hasNext() {
            return index < size;
        }

        public VectorEntry next() {
            entry.update(index);

            index++;

            return entry;
        }

        public void remove() {
            entry.set(0);
        }

    }

    /**
     * Vector entry backed by the vector. May be reused for higher performance
     */
    private class RefVectorEntry implements VectorEntry {

        private int index;

        /**
         * Updates the entry
         */
        public void update(int index) {
            this.index = index;
        }

        public int index() {
            return index;
        }

        public double get() {
            return AbstractVector.this.get(index);
        }

        public void set(double value) {
            AbstractVector.this.set(index, value);
        }

    }

}
