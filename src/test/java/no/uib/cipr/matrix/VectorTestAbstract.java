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

import no.uib.cipr.matrix.Vector.Norm;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Test of vectors
 */
public abstract class VectorTestAbstract {

    /**
     * Max vector size
     */
    protected int max = 100;

    /**
     * Vector to test
     */
    protected Vector x;

    /**
     * Data of the vector
     */
    protected double[] xd;

    /**
     * Vectors of the same size as x, dense and non-dense
     */
    protected Vector y, yDense, z, zDense;

    /**
     * Data of those vectors
     */
    protected double[] yd, zd;

    /**
     * Tolerance for floating-point comparisons
     */
    protected double tol = 1e-5;

    @Before
    public void setUp() throws Exception {
        createPrimary();
        createAuxillerary();
    }

    protected abstract void createPrimary() throws Exception;

    protected void createAuxillerary() throws Exception {
        yDense = Matrices.random(x.size());
        y = Matrices.synchronizedVector(yDense);
        zDense = Matrices.random(x.size());
        z = Matrices.synchronizedVector(zDense);
        yd = Matrices.getArray(y);
        zd = Matrices.getArray(z);
    }

    @After
    public void tearDown() throws Exception {
        x = null;
        xd = null;
    }

    @Test
    public void testSize() {
        assertEquals(xd.length, x.size());
    }

    @Test
    public void testCopy() {
        Vector y = x.copy();
        assertVectorEquals(xd, y);
        y.scale(Math.random());
        assertVectorEquals(xd, x);
    }

    @Test
    public void testZero() {
        assertVectorEquals(scale(0), x.zero());
    }

    @Test
    public void testScale() {
        double alpha = Math.random();
        assertVectorEquals(scale(alpha), x.scale(alpha));
    }

    @Test
    public void testScaleZero() {
        assertVectorEquals(scale(0), x.scale(0));
    }

    @Test
    public void testScaleOne() {
        assertVectorEquals(scale(1), x.scale(1));
    }

    protected double[] scale(double alpha) {
        for (int i = 0; i < xd.length; ++i)
            xd[i] *= alpha;
        return xd;
    }

    /*
     * Test for Vector set(Vector)
     */
    @Test
    public void testSetVectorDense() {
        assertVectorEquals(set(1, yd, 0, zd), x.set(yDense));
    }

    /*
     * Test for Vector set(double, Vector)
     */
    @Test
    public void testSetdoubleVectorDense() {
        double alpha = Math.random();
        assertVectorEquals(set(alpha, yd, 0, zd), x.set(alpha, yDense));
    }

    /*
     * Test for Vector set(Vector)
     */
    @Test
    public void testSetVector() {
        assertVectorEquals(set(1, yd, 0, zd), x.set(y));
    }

    /*
     * Test for Vector set(double, Vector)
     */
    @Test
    public void testSetDoubleVector() {
        double alpha = Math.random();
        assertVectorEquals(set(alpha, yd, 0, zd), x.set(alpha, y));
    }

    protected double[] set(double alpha, double[] yd, double beta, double[] zd) {
        for (int i = 0; i < xd.length; ++i)
            xd[i] = alpha * yd[i] + beta * zd[i];
        return xd;
    }

    /*
     * Test for Vector add(double, Vector)
     */
    @Test
    public void testAddDoubleVectorDense() {
        double alpha = Math.random();
        assertVectorEquals(add(alpha, yd), x.add(alpha, yDense));
    }

    /*
     * Test for Vector add(Vector)
     */
    @Test
    public void testAddVectorDense() {
        assertVectorEquals(add(1, yd), x.add(yDense));
    }

    /*
     * Test for Vector add(double, Vector)
     */
    @Test
    public void testAddDoubleVector() {
        double alpha = Math.random();
        assertVectorEquals(add(alpha, yd), x.add(alpha, y));
    }

    /*
     * Test for Vector add(Vector)
     */
    @Test
    public void testAddVector() {
        assertVectorEquals(add(1, yd), x.add(y));
    }

    protected double[] add(double alpha, double[] yd) {
        for (int i = 0; i < xd.length; ++i)
            xd[i] += alpha * yd[i];
        return xd;
    }

    @Test
    public void testDotDense() {
        assertEquals(dot(yd), x.dot(yDense), tol);
        assertVectorEquals(xd, x);
        assertVectorEquals(yd, yDense);
    }

    @Test
    public void testDot() {
        assertEquals(dot(yd), x.dot(y), tol);
        assertVectorEquals(xd, x);
        assertVectorEquals(yd, y);
    }

    protected double dot(double[] yd) {
        double dot = 0;
        for (int i = 0; i < yd.length; ++i)
            dot += xd[i] * yd[i];
        return dot;
    }

    @Test
    public void testCardinality() {
        assertEquals(cardinality(), Matrices.cardinality(x));
    }

    protected int cardinality() {
        int nz = 0;
        for (double d : xd)
            if (d != 0)
                nz++;
        return nz;
    }

    @Test
    public void testNorm1() {
        assertEquals(norm1(), x.norm(Norm.One), tol);
    }

    @Test
    public void testNorm2() {
        assertEquals(norm2(), x.norm(Norm.Two), tol);
    }

    @Test
    public void testNormInf() {
        assertEquals(normInf(), x.norm(Norm.Infinity), tol);
    }

    protected double norm1() {
        double norm = 0;
        for (double d : xd)
            norm += Math.abs(d);
        return norm;
    }

    protected double norm2() {
        double norm = 0;
        for (double d : xd)
            norm += d * d;
        return Math.sqrt(norm);
    }

    protected double normInf() {
        double norm = 0;
        for (double d : xd)
            norm = Math.max(Math.abs(d), norm);
        return norm;
    }

    @Test
    public void testIteratorSetGet() {
        double alpha = Math.random();
        double[] data = new double[x.size()];
        for (VectorEntry e : x) {
            data[e.index()] = e.get();
            e.set(alpha * e.get());
            e.set(e.get() / alpha);
        }
        assertVectorEquals(xd, x);
        assertVectorEquals(xd, data);
    }

    @Test
    public void testIteratorGet() {
        double[] data = new double[x.size()];
        for (VectorEntry e : x)
            data[e.index()] = e.get();
        assertVectorEquals(xd, x);
        assertVectorEquals(xd, data);
    }

    @Test
    public void testIteratorSet() {
        double alpha = Math.random();
        for (VectorEntry e : x)
            e.set(alpha * e.get());
        assertVectorEquals(scale(alpha), x);
    }

    protected void assertVectorEquals(double[] xd, Vector x) {
        assertEquals(xd.length, x.size());
        for (int i = 0; i < xd.length; ++i)
            assertEquals(xd[i], x.get(i), tol);
    }

    protected void assertVectorEquals(double[] xd, double[] yd) {
        for (int i = 0; i < xd.length; ++i)
            assertEquals(xd[i], yd[i], tol);
    }

}
