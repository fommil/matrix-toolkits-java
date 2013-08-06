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

import junit.framework.TestCase;
import no.uib.cipr.matrix.Vector.Norm;

/**
 * Test of vectors
 */
public abstract class VectorTestAbstract extends TestCase {

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

    public VectorTestAbstract(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
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

    @Override
    protected void tearDown() throws Exception {
        x = null;
        xd = null;
    }

    public void testSize() {
        assertEquals(xd.length, x.size());
    }

    public void testCopy() {
        Vector y = x.copy();
        assertEquals(xd, y);
        y.scale(Math.random());
        assertEquals(xd, x);
    }

    public void testZero() {
        assertEquals(scale(0), x.zero());
    }

    public void testScale() {
        double alpha = Math.random();
        assertEquals(scale(alpha), x.scale(alpha));
    }

    public void testScaleZero() {
        assertEquals(scale(0), x.scale(0));
    }

    public void testScaleOne() {
        assertEquals(scale(1), x.scale(1));
    }

    protected double[] scale(double alpha) {
        for (int i = 0; i < xd.length; ++i)
            xd[i] *= alpha;
        return xd;
    }

    /*
     * Test for Vector set(Vector)
     */
    public void testSetVectorDense() {
        assertEquals(set(1, yd, 0, zd), x.set(yDense));
    }

    /*
     * Test for Vector set(double, Vector)
     */
    public void testSetdoubleVectorDense() {
        double alpha = Math.random();
        assertEquals(set(alpha, yd, 0, zd), x.set(alpha, yDense));
    }

    /*
     * Test for Vector set(Vector)
     */
    public void testSetVector() {
        assertEquals(set(1, yd, 0, zd), x.set(y));
    }

    /*
     * Test for Vector set(double, Vector)
     */
    public void testSetDoubleVector() {
        double alpha = Math.random();
        assertEquals(set(alpha, yd, 0, zd), x.set(alpha, y));
    }

    protected double[] set(double alpha, double[] yd, double beta, double[] zd) {
        for (int i = 0; i < xd.length; ++i)
            xd[i] = alpha * yd[i] + beta * zd[i];
        return xd;
    }

    /*
     * Test for Vector add(double, Vector)
     */
    public void testAddDoubleVectorDense() {
        double alpha = Math.random();
        assertEquals(add(alpha, yd), x.add(alpha, yDense));
    }

    /*
     * Test for Vector add(Vector)
     */
    public void testAddVectorDense() {
        assertEquals(add(1, yd), x.add(yDense));
    }

    /*
     * Test for Vector add(double, Vector)
     */
    public void testAddDoubleVector() {
        double alpha = Math.random();
        assertEquals(add(alpha, yd), x.add(alpha, y));
    }

    /*
     * Test for Vector add(Vector)
     */
    public void testAddVector() {
        assertEquals(add(1, yd), x.add(y));
    }

    protected double[] add(double alpha, double[] yd) {
        for (int i = 0; i < xd.length; ++i)
            xd[i] += alpha * yd[i];
        return xd;
    }

    public void testDotDense() {
        assertEquals(dot(yd), x.dot(yDense), tol);
        assertEquals(xd, x);
        assertEquals(yd, yDense);
    }

    public void testDot() {
        assertEquals(dot(yd), x.dot(y), tol);
        assertEquals(xd, x);
        assertEquals(yd, y);
    }

    protected double dot(double[] yd) {
        double dot = 0;
        for (int i = 0; i < yd.length; ++i)
            dot += xd[i] * yd[i];
        return dot;
    }

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

    public void testNorm1() {
        assertEquals(norm1(), x.norm(Norm.One), tol);
    }

    public void testNorm2() {
        assertEquals(norm2(), x.norm(Norm.Two), tol);
    }

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

    public void testIteratorSetGet() {
        double alpha = Math.random();
        double[] data = new double[x.size()];
        for (VectorEntry e : x) {
            data[e.index()] = e.get();
            e.set(alpha * e.get());
            e.set(e.get() / alpha);
        }
        assertEquals(xd, x);
        assertEquals(xd, data);
    }

    public void testIteratorGet() {
        double[] data = new double[x.size()];
        for (VectorEntry e : x)
            data[e.index()] = e.get();
        assertEquals(xd, x);
        assertEquals(xd, data);
    }

    public void testIteratorSet() {
        double alpha = Math.random();
        for (VectorEntry e : x)
            e.set(alpha * e.get());
        assertEquals(scale(alpha), x);
    }

    protected void assertEquals(double[] xd, Vector x) {
        assertEquals(xd.length, x.size());
        for (int i = 0; i < xd.length; ++i)
            assertEquals(xd[i], x.get(i), tol);
    }

    protected void assertEquals(double[] xd, double[] yd) {
        for (int i = 0; i < xd.length; ++i)
            assertEquals(xd[i], yd[i], tol);
    }

}
