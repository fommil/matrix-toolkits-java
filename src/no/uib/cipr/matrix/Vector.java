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

/**
 * Basic vector interface. It holds <code>double</code>s in an array, and is
 * used alongside <code>Matrix</code> in numerical computations. Implementing
 * classes decides on the actual storage.
 * 
 * <h4>Basic operations</h4>
 * <p>
 * Use <code>size</code> to get the vector size. <code>get(int)</code> gets
 * an element, and there are corresponding <code>set(int,double)</code> and
 * <code>add(int,double)</code> methods as well. Note that vector indices are
 * zero-based (typical for Java and C). This means that they range from 0 to
 * <code>size-1</code>. It is legal to have <code>size</code> equal zero.
 * </p>
 * <p>
 * Other basic operations are <code>zero</code> which zeros all the entries of
 * the vector, which can be cheaper than either zeroing the vector manually, or
 * creating a new vector, and the operation <code>copy</code> which creates a
 * deep copy of the vector. This copy has separate storage, but starts with the
 * same contents as the current vector.
 * </p>
 * <h4>Iterators</h4>
 * <p>
 * The vector interface extends <code>Iterable</code>, and the iterator
 * returns a <code>VectorEntry</code> which contains current index and entry
 * value. Note that the iterator may skip non-zero entries. Using an iterator,
 * many simple and efficient algorithms can be created. The iterator also
 * permits changing values in the vector, however only non-zero entries can be
 * changed.
 * </p>
 * <h4>Basic linear algebra</h4>
 * <p>
 * A selection of basic linear algebra operations are available. To ensure high
 * efficiency, little or no internal memory allocation is done, and the user is
 * required to supply the output arguments.
 * </p>
 * <p>
 * The operations available include:
 * </p>
 * <dl>
 * <dt><i>Additions</i></dt>
 * <dd>Vectors can be added to each other, even if their underlying vector
 * structures are incompatible</dd>
 * <dt><i>Scaling</i></dt>
 * <dd>Scalar multiplication (scaling) of a whole vector</dd>
 * <dt><i>Norms</i></dt>
 * <dd>Both innerproducts and norms can be computed. Several common norms are
 * supported</dd>
 * </dl>
 */
public interface Vector extends Iterable<VectorEntry>, Serializable {

    /**
     * Size of the vector
     */
    int size();

    /**
     * <code>x(index) = value</code>
     */
    void set(int index, double value);

    /**
     * <code>x(index) += value</code>
     */
    void add(int index, double value);

    /**
     * Returns <code>x(index)</code>
     */
    double get(int index);

    /**
     * Creates a deep copy of the vector
     */
    Vector copy();

    /**
     * Zeros all the entries in the vector, while preserving any underlying
     * structure
     */
    Vector zero();

    /**
     * <code>x=alpha*x</code>
     * 
     * @return x
     */
    Vector scale(double alpha);

    /**
     * <code>x=y</code>
     * 
     * @return x
     */
    Vector set(Vector y);

    /**
     * <code>x=alpha*y</code>
     * 
     * @return x
     */
    Vector set(double alpha, Vector y);

    /**
     * <code>x = y + x</code>
     * 
     * @return x
     */
    Vector add(Vector y);

    /**
     * <code>x = alpha*y + x</code>
     * 
     * @return x
     */
    Vector add(double alpha, Vector y);

    /**
     * <code>x<sup>T</sup>*y</code>
     */
    double dot(Vector y);

    /**
     * Computes the given norm of the vector
     * 
     * @param type
     *            The type of norm to compute
     */
    double norm(Norm type);

    /**
     * Supported vector-norms. The difference between the two 2-norms is that
     * one is fast, but can overflow, while the robust version is overflow
     * resistant, but slower.
     */
    enum Norm {

        /**
         * Sum of the absolute values of the entries
         */
        One,

        /**
         * The root of sum of squares
         */
        Two,

        /**
         * As the 2 norm may overflow, an overflow resistant version is also
         * available. Note that it may be slower.
         */
        TwoRobust,

        /**
         * Largest entry in absolute value
         */
        Infinity

    }

}
