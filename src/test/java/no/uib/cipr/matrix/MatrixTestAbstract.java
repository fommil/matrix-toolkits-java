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

/**
 * Tests a matrix
 */
public abstract class MatrixTestAbstract extends TestCase {

    /**
     * Matrix to test
     */
    protected Matrix A;

    /**
     * Jagged array version of A
     */
    protected double[][] Ad;

    /**
     * Matrix of the same size as A, dense and non-dense
     */
    protected Matrix Bdense, B;

    /**
     * Contents of B
     */
    protected double[][] Bd;

    /**
     * Non-dense vectors with size equal the number of rows in A
     */
    protected Vector xR, yR;

    /**
     * Non-dense vectors with size equal the number of columns in A
     */
    protected Vector xC, yC;

    /**
     * Dense vectors with size equal the number of rows in A
     */
    protected Vector xDenseR, yDenseR;

    /**
     * Dense vectors with size equal the number of columns in A
     */
    protected Vector xDenseC, yDenseC;

    /**
     * Contents of the vectors
     */
    protected double[] xdR, ydR, xdC, ydC;

    /**
     * Tolerance for floating-point comparisons
     */
    protected double tol = 1e-4;

    /**
     * Maximum matrix size, to avoid too slow tests
     */
    protected int max = 100;

    /**
     * Constructor for MatrixTestAbstract
     */
    public MatrixTestAbstract(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
        createPrimary();
        createAuxillerary();
    }

    protected abstract void createPrimary() throws Exception;

    @Override
    protected void tearDown() throws Exception {
        A = B = Bdense = null;
        Ad = Bd = null;
        xC = xDenseC = xDenseR = xR = yC = yDenseC = yDenseR = yR = null;
        xdC = xdR = ydC = ydR = null;
    }

    /**
     * Called after setUp() to create additional datastructures
     */
    protected void createAuxillerary() {
        Bdense = Matrices.random(A.numRows(), A.numColumns());
        B = Matrices.synchronizedMatrix(Bdense.copy());
        Bd = Matrices.getArray(B);

        xDenseC = Matrices.random(A.numColumns());
        yDenseC = Matrices.random(A.numColumns());

        xDenseR = Matrices.random(A.numRows());
        yDenseR = Matrices.random(A.numRows());

        xC = Matrices.synchronizedVector(xDenseC);
        yC = Matrices.synchronizedVector(yDenseC);

        xR = Matrices.synchronizedVector(xDenseR);
        yR = Matrices.synchronizedVector(yDenseR);

        xdC = Matrices.getArray(xC);
        ydC = Matrices.getArray(yC);

        xdR = Matrices.getArray(xR);
        ydR = Matrices.getArray(yR);
    }

    public void testMatrixRank2Dense() {
        if (A.isSquare()) {
            int n = Utilities.getInt(1, max);
            Matrix B = Matrices.random(A.numRows(), n), C = Matrices.random(A
                    .numRows(), n);
            double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
            double alpha = Math.random();

            A = A.rank2(alpha, B, C);
            rank2(Ad, alpha, Bd, Cd);

            assertEquals(Ad, A);
            assertEquals(Bd, B);
            assertEquals(Cd, C);
        }
    }

    public void testMatrixRank2() {
        if (A.isSquare()) {
            int n = Utilities.getInt(1, max);
            Matrix B = Matrices.synchronizedMatrix(Matrices.random(A.numRows(),
                    n)), C = Matrices.synchronizedMatrix(Matrices.random(A
                    .numRows(), n));
            double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
            double alpha = Math.random();

            A = A.rank2(alpha, B, C);
            rank2(Ad, alpha, Bd, Cd);

            assertEquals(Ad, A);
            assertEquals(Bd, B);
            assertEquals(Cd, C);
        }
    }

    public void testMatrixTransRank2Dense() {
        if (A.isSquare()) {
            int n = Utilities.getInt(1, max);
            Matrix B = Matrices.random(n, A.numColumns()), C = Matrices.random(
                    n, A.numColumns());
            double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
            double alpha = Math.random();

            A = A.transRank2(alpha, B, C);
            transRank2(Ad, alpha, Bd, Cd);

            assertEquals(Ad, A);
            assertEquals(Bd, B);
            assertEquals(Cd, C);
        }
    }

    public void testMatrixTransRank2() {
        if (A.isSquare()) {
            int n = Utilities.getInt(1, max);
            Matrix B = Matrices.synchronizedMatrix(Matrices.random(n, A
                    .numColumns())), C = Matrices.synchronizedMatrix(Matrices
                    .random(n, A.numColumns()));
            double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
            double alpha = Math.random();

            A = A.transRank2(alpha, B, C);
            transRank2(Ad, alpha, Bd, Cd);

            assertEquals(Ad, A);
            assertEquals(Bd, B);
            assertEquals(Cd, C);
        }
    }

    public void testMatrixRank1Dense() {
        if (A.isSquare()) {
            Matrix C = Matrices.random(A.numRows(), A.numColumns());
            double[][] Cd = Matrices.getArray(C);
            double alpha = Math.random();

            A = A.rank1(alpha, C);
            rank1(Ad, alpha, Cd);

            assertEquals(Ad, A);
            assertEquals(Cd, C);
        }
    }

    public void testMatrixRank1() {
        if (A.isSquare()) {
            Matrix C = Matrices.synchronizedMatrix(Matrices.random(A.numRows(),
                    A.numColumns()));
            double[][] Cd = Matrices.getArray(C);
            double alpha = Math.random();

            A = A.rank1(alpha, C);
            rank1(Ad, alpha, Cd);

            assertEquals(Ad, A);
            assertEquals(Cd, C);
        }
    }

    public void testMatrixTransRank1Dense() {
        if (A.isSquare()) {
            Matrix C = Matrices.random(A.numRows(), A.numColumns());
            double[][] Cd = Matrices.getArray(C);
            double alpha = Math.random();

            A = A.transRank1(alpha, C);
            transRank1(Ad, alpha, Cd);

            assertEquals(Ad, A);
            assertEquals(Cd, C);
        }
    }

    public void testMatrixTransRank1() {
        if (A.isSquare()) {
            Matrix C = Matrices.synchronizedMatrix(Matrices.random(A.numRows(),
                    A.numColumns()));
            double[][] Cd = Matrices.getArray(C);
            double alpha = Math.random();

            A = A.transRank1(alpha, C);
            transRank1(Ad, alpha, Cd);

            assertEquals(Ad, A);
            assertEquals(Cd, C);
        }
    }

    public void testMatrixMultDense() {
        int m = A.numRows(), k = A.numColumns(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.random(k, n), C = Matrices.random(m, n);
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.mult(alpha, B, C);
        Cd = mult(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixMult() {
        int m = A.numRows(), k = A.numColumns(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.synchronizedMatrix(Matrices.random(k, n)), C = Matrices
                .synchronizedMatrix(Matrices.random(m, n));
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.mult(alpha, B, C);
        Cd = mult(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransAmultDense() {
        int m = A.numColumns(), k = A.numRows(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.random(k, n), C = Matrices.random(m, n);
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transAmult(alpha, B, C);
        Cd = transAmult(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransAmult() {
        int m = A.numColumns(), k = A.numRows(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.synchronizedMatrix(Matrices.random(k, n)), C = Matrices
                .synchronizedMatrix(Matrices.random(m, n));
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transAmult(alpha, B, C);
        Cd = transAmult(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransABmultDense() {
        int m = A.numColumns(), k = A.numRows(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.random(n, k), C = Matrices.random(m, n);
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transABmult(alpha, B, C);
        Cd = transABmult(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransABmult() {
        int m = A.numColumns(), k = A.numRows(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.synchronizedMatrix(Matrices.random(n, k)), C = Matrices
                .synchronizedMatrix(Matrices.random(m, n));
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transABmult(alpha, B, C);
        Cd = transABmult(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransBmultDense() {
        int m = A.numRows(), k = A.numColumns(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.random(n, k), C = Matrices.random(m, n);
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transBmult(alpha, B, C);
        Cd = transBmult(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransBmult() {
        int m = A.numRows(), k = A.numColumns(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.synchronizedMatrix(Matrices.random(n, k)), C = Matrices
                .synchronizedMatrix(Matrices.random(m, n));
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transBmult(alpha, B, C);
        Cd = transBmult(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixMultAddDense() {
        int m = A.numRows(), k = A.numColumns(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.random(k, n), C = Matrices.random(m, n);
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.multAdd(alpha, B, C);
        Cd = multAdd(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixMultAdd() {
        int m = A.numRows(), k = A.numColumns(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.synchronizedMatrix(Matrices.random(k, n)), C = Matrices
                .synchronizedMatrix(Matrices.random(m, n));
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.multAdd(alpha, B, C);
        Cd = multAdd(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransAmultAddDense() {
        int m = A.numColumns(), k = A.numRows(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.random(k, n), C = Matrices.random(m, n);
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transAmultAdd(alpha, B, C);
        Cd = transAmultAdd(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransAmultAdd() {
        int m = A.numColumns(), k = A.numRows(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.synchronizedMatrix(Matrices.random(k, n)), C = Matrices
                .synchronizedMatrix(Matrices.random(m, n));
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transAmultAdd(alpha, B, C);
        Cd = transAmultAdd(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransABmultAddDense() {
        int m = A.numColumns(), k = A.numRows(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.random(n, k), C = Matrices.random(m, n);
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transABmultAdd(alpha, B, C);
        Cd = transABmultAdd(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransABmultAdd() {
        int m = A.numColumns(), k = A.numRows(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.synchronizedMatrix(Matrices.random(n, k)), C = Matrices
                .synchronizedMatrix(Matrices.random(m, n));
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transABmultAdd(alpha, B, C);
        Cd = transABmultAdd(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransBmultAddDense() {
        int m = A.numRows(), k = A.numColumns(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.random(n, k), C = Matrices.random(m, n);
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transBmultAdd(alpha, B, C);
        Cd = transBmultAdd(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    public void testMatrixTransBmultAdd() {
        int m = A.numRows(), k = A.numColumns(), n = Utilities.getInt(1, max);
        Matrix B = Matrices.synchronizedMatrix(Matrices.random(n, k)), C = Matrices
                .synchronizedMatrix(Matrices.random(m, n));
        double[][] Bd = Matrices.getArray(B), Cd = Matrices.getArray(C);
        double alpha = Math.random();

        C = A.transBmultAdd(alpha, B, C);
        Cd = transBmultAdd(Ad, alpha, Bd, Cd);

        assertEquals(Ad, A);
        assertEquals(Bd, B);
        assertEquals(Cd, C);
    }

    protected double[][] rank2(double[][] Ad, double alpha, double[][] Bd,
            double[][] Cd) {
        return transBmultAdd(Bd, alpha, Cd, transBmultAdd(Cd, alpha, Bd, Ad));
    }

    protected double[][] transRank2(double[][] Ad, double alpha, double[][] Bd,
            double[][] Cd) {
        return transAmultAdd(Bd, alpha, Cd, transAmultAdd(Cd, alpha, Bd, Ad));
    }

    protected double[][] rank1(double[][] Ad, double alpha, double[][] Cd) {
        return transBmultAdd(Cd, alpha, Cd, Ad);
    }

    protected double[][] transRank1(double[][] Ad, double alpha, double[][] Cd) {
        return transAmultAdd(Cd, alpha, Cd, Ad);
    }

    public void testVectorRank2Dense() {
        if (A.isSquare()) {
            double alpha = Math.random();
            assertEquals(rank2(alpha, xdR, ydR), A.rank2(alpha, xDenseR,
                    yDenseR));
        }
    }

    public void testVectorRank2() {
        if (A.isSquare()) {
            double alpha = Math.random();
            assertEquals(rank2(alpha, xdR, ydR), A.rank2(alpha, xR, yR));
        }
    }

    public void testVectorRank1Dense() {
        if (A.isSquare()) {
            double alpha = Math.random();
            assertEquals(rank1(alpha, xdR, ydR), A.rank1(alpha, xDenseR,
                    yDenseR));
        }
    }

    public void testVectorRank1() {
        if (A.isSquare()) {
            double alpha = Math.random();
            assertEquals(rank1(alpha, xdR, ydR), A.rank1(alpha, xR, yR));
        }
    }

    protected double[][] rank2(double alpha, double[] xd, double[] yd) {
        rank1(alpha, xd, yd);
        rank1(alpha, yd, xd);
        return Ad;
    }

    protected double[][] rank1(double alpha, double[] xd, double[] yd) {
        for (int i = 0; i < xd.length; ++i)
            for (int j = 0; j < yd.length; ++j)
                Ad[i][j] += alpha * xd[i] * yd[j];
        return Ad;
    }

    public void testVectorTransMultAddDense() {
        double alpha = Math.random();
        assertEquals(transMultAdd(alpha, xdR, ydC), A.transMultAdd(alpha,
                xDenseR, yDenseC));
        assertEquals(Ad, A);
        assertEquals(xdR, xDenseR);
        assertEquals(ydC, yDenseC);
    }

    public void testVectorTransMultAdd() {
        double alpha = Math.random();
        assertEquals(transMultAdd(alpha, xdR, ydC), A.transMultAdd(alpha, xR,
                yC));
        assertEquals(Ad, A);
        assertEquals(xdR, xR);
        assertEquals(ydC, yC);
    }

    protected double[] transMultAdd(double alpha, double[] xd, double[] yd) {
        int rows = Ad.length, cols = 0;
        if (rows > 0)
            cols = Ad[0].length;

        for (int j = 0; j < cols; ++j) {
            double dot = 0;
            for (int i = 0; i < rows; ++i)
                dot += Ad[i][j] * xd[i];
            yd[j] += alpha * dot;
        }

        return yd;
    }

    public void testVectorMultDense() {
        double alpha = Math.random();
        assertEquals(mult(alpha, xdC, ydR), A.mult(alpha, xDenseC, yDenseR));
        assertEquals(Ad, A);
        assertEquals(xdC, xDenseC);
        assertEquals(ydR, yDenseR);
    }

    public void testVectorMult() {
        double alpha = Math.random();
        assertEquals(mult(alpha, xdC, ydR), A.mult(alpha, xC, yR));
        assertEquals(Ad, A);
        assertEquals(xdC, xC);
        assertEquals(ydR, yR);
    }

    protected double[] mult(double alpha, double[] xd, double[] yd) {
        for (int i = 0; i < Ad.length; ++i) {
            double dot = 0;
            for (int j = 0; j < Ad[i].length; ++j)
                dot += Ad[i][j] * xd[j];
            yd[i] = alpha * dot;
        }
        return yd;
    }

    public void testVectorMultAddDense() {
        double alpha = Math.random();
        assertEquals(multAdd(Ad, alpha, xdC, ydR), A.multAdd(alpha, xDenseC,
                yDenseR));
        assertEquals(Ad, A);
        assertEquals(xdC, xDenseC);
        assertEquals(ydR, yDenseR);
    }

    public void testVectorMultAdd() {
        double alpha = Math.random();
        assertEquals(multAdd(Ad, alpha, xdC, ydR), A.multAdd(alpha, xC, yR));
        assertEquals(Ad, A);
        assertEquals(xdC, xC);
        assertEquals(ydR, yR);
    }

    protected double[] multAdd(double[][] Ad, double alpha, double[] xd,
            double[] yd) {
        for (int i = 0; i < Ad.length; ++i) {
            double dot = 0;
            for (int j = 0; j < Ad[i].length; ++j)
                dot += Ad[i][j] * xd[j];
            yd[i] += alpha * dot;
        }

        return yd;
    }

    protected double[][] mult(double[][] Ad, double alpha, double[][] Bd,
            double[][] Cd) {
        int m = Cd.length, n = 0, k = Bd.length;
        if (k > 0)
            n = Bd[0].length;

        Utilities.zero(Cd);

        for (int j = 0; j < n; ++j)
            for (int l = 0; l < k; ++l)
                for (int i = 0; i < m; ++i)
                    Cd[i][j] += alpha * Ad[i][l] * Bd[l][j];

        return Cd;
    }

    protected double[][] transAmult(double[][] Ad, double alpha, double[][] Bd,
            double[][] Cd) {
        int m = Cd.length, n = 0, k = Bd.length;
        if (k > 0)
            n = Bd[0].length;

        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i) {
                double temp = 0;
                for (int l = 0; l < k; ++l)
                    temp += Ad[l][i] * Bd[l][j];
                Cd[i][j] = alpha * temp;
            }

        return Cd;
    }

    protected double[][] transBmult(double[][] Ad, double alpha, double[][] Bd,
            double[][] Cd) {
        int m = Cd.length, n = Bd.length, k = 0;
        if (n > 0)
            k = Bd[0].length;

        Utilities.zero(Cd);

        for (int j = 0; j < n; ++j) {
            for (int l = 0; l < k; ++l)
                for (int i = 0; i < m; ++i)
                    Cd[i][j] += alpha * Ad[i][l] * Bd[j][l];
        }

        return Cd;
    }

    protected double[][] transABmult(double[][] Ad, double alpha,
            double[][] Bd, double[][] Cd) {
        int m = Cd.length, n = Bd.length, k = 0;
        if (n > 0)
            k = Bd[0].length;

        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i) {
                double temp = 0;
                for (int l = 0; l < k; ++l)
                    temp += Ad[l][i] * Bd[j][l];
                Cd[i][j] = alpha * temp;
            }

        return Cd;
    }

    protected double[][] multAdd(double[][] Ad, double alpha, double[][] Bd,
            double[][] Cd) {
        int m = Cd.length, n = 0, k = Bd.length;
        if (k > 0)
            n = Bd[0].length;

        for (int j = 0; j < n; ++j)
            for (int l = 0; l < k; ++l)
                for (int i = 0; i < m; ++i)
                    Cd[i][j] += alpha * Ad[i][l] * Bd[l][j];

        return Cd;
    }

    protected double[][] transAmultAdd(double[][] Ad, double alpha,
            double[][] Bd, double[][] Cd) {
        int m = Cd.length, n = 0, k = Bd.length;
        if (k > 0)
            n = Bd[0].length;

        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i) {
                double temp = 0;
                for (int l = 0; l < k; ++l)
                    temp += Ad[l][i] * Bd[l][j];
                Cd[i][j] += alpha * temp;
            }

        return Cd;
    }

    protected double[][] transBmultAdd(double[][] Ad, double alpha,
            double[][] Bd, double[][] Cd) {
        int m = Cd.length, n = Bd.length, k = 0;
        if (n > 0)
            k = Bd[0].length;

        for (int j = 0; j < n; ++j)
            for (int l = 0; l < k; ++l)
                for (int i = 0; i < m; ++i)
                    Cd[i][j] += alpha * Ad[i][l] * Bd[j][l];

        return Cd;
    }

    protected double[][] transABmultAdd(double[][] Ad, double alpha,
            double[][] Bd, double[][] Cd) {
        int m = Cd.length, n = Bd.length, k = 0;
        if (n > 0)
            k = Bd[0].length;

        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i) {
                double temp = 0;
                for (int l = 0; l < k; ++l)
                    temp += Ad[l][i] * Bd[j][l];
                Cd[i][j] += alpha * temp;
            }

        return Cd;
    }

    /**
     * Tests <code>A = A + alpha*B</code>
     */
    public void testRandomMatrixAdd() {
        double alpha = Math.random();
        A = A.add(alpha, B);
        add(Ad, alpha, Bd);
        assertEquals(Ad, A);
        assertEquals(Bd, B);
    }

    /**
     * Tests <code>A = A + B</code>
     */
    public void testMatrixAdd() {
        A = A.add(B);
        add(Ad, 1, Bd);
        assertEquals(Ad, A);
        assertEquals(Bd, B);
    }

    /**
     * Tests <code>A = A + 1*B</code>
     */
    public void testOneMatrixAdd() {
        A = A.add(1, B);
        add(Ad, 1, Bd);
        assertEquals(Ad, A);
        assertEquals(Bd, B);
    }

    /**
     * Tests <code>A = A + 0*B</code>
     */
    public void testZeroMatrixAdd() {
        A = A.add(0, B);
        add(Ad, 0, Bd);
        assertEquals(Ad, A);
        assertEquals(Bd, B);
    }

    /**
     * Tests <code>A = alpha*B</code>
     */
    public void testRandomMatrixSet() {
        double alpha = Math.random();
        A = A.set(alpha, B);
        set(Ad, alpha, Bd);
        assertEquals(Ad, A);
        assertEquals(Bd, B);
    }

    /**
     * Tests <code>A = B</code>
     */
    public void testMatrixSet() {
        A = A.set(B);
        set(Ad, 1, Bd);
        assertEquals(Ad, A);
        assertEquals(Bd, B);
    }

    /**
     * Tests <code>A = 1*B</code>
     */
    public void testOneMatrixSet() {
        A = A.set(1, B);
        set(Ad, 1, Bd);
        assertEquals(Ad, A);
        assertEquals(Bd, B);
    }

    /**
     * Tests <code>A = 0*B</code>
     */
    public void testZeroMatrixSet() {
        A = A.set(0, B);
        set(Ad, 0, Bd);
        assertEquals(Ad, A);
        assertEquals(Bd, B);
    }

    /**
     * Checks transpose
     */
    public void testTranspose() {
        Matrix At = Matrices.random(A.numColumns(), A.numRows());
        assertEquals(transpose(), A.transpose(At));
    }

    protected void set(double[][] A, double alpha, double[][] B) {
        for (int i = 0; i < A.length; ++i)
            for (int j = 0; j < A[i].length; ++j)
                A[i][j] = alpha * B[i][j];
    }

    protected void add(double[][] A, double alpha, double[][] B) {
        for (int i = 0; i < A.length; ++i)
            for (int j = 0; j < A[i].length; ++j)
                A[i][j] += alpha * B[i][j];
    }

    protected double[][] transpose() {
        if (Ad.length == 0)
            return new double[0][0];

        double[][] Adt = new double[Ad[0].length][Ad.length];
        for (int i = 0; i < Ad.length; ++i)
            for (int j = 0; j < Ad[i].length; ++j)
                Adt[j][i] = Ad[i][j];
        return Adt;
    }

    /**
     * Test of direct matrix solver
     */
    public void testMatrixSolve() {
        while (true) {
            try {
                Matrix B = Matrices.random(A.numRows(), A.numColumns());
                Matrix X = Matrices.random(A.numRows(), A.numColumns());
                X = A.solve(B, X);

                Matrix Y = A.multAdd(X, X.copy().set(-1, B));
                assertEquals(0, Y.norm(Matrix.Norm.Frobenius), tol);
                assertEquals(Ad, A);
                return;
            } catch (MatrixSingularException e) {
                Utilities.addDiagonal(A, Ad, 1);
            } catch (MatrixNotSPDException e) {
                Utilities.addDiagonal(A, Ad, 1);
            }
        }
    }

    /**
     * Test of direct transpose matrix solver
     */
    public void testTransMatrixSolve() {
        while (true) {
            try {
                Matrix B = Matrices.random(A.numRows(), A.numColumns());
                Matrix X = Matrices.random(A.numRows(), A.numColumns());
                X = A.transSolve(B, X);

                Matrix Y = A.transAmultAdd(X, X.copy().set(-1, B));
                assertEquals(0, Y.norm(Matrix.Norm.Frobenius), tol);
                assertEquals(Ad, A);
                return;
            } catch (MatrixSingularException e) {
                Utilities.addDiagonal(A, Ad, 1);
            } catch (MatrixNotSPDException e) {
                Utilities.addDiagonal(A, Ad, 1);
            }
        }
    }

    /**
     * Test of direct vector solver
     */
    public void testVectorSolve() {
        while (true) {
            try {
                Vector b = Matrices.random(A.numRows());
                Vector x = Matrices.random(A.numRows());
                x = A.solve(b, x);

                Vector y = A.multAdd(-1, x, x.copy().set(b));
                assertEquals(0, y.norm(Vector.Norm.Two), tol);
                assertEquals(Ad, A);
                return;
            } catch (MatrixSingularException e) {
                Utilities.addDiagonal(A, Ad, 1);
            } catch (MatrixNotSPDException e) {
                Utilities.addDiagonal(A, Ad, 1);
            }
        }
    }

    /**
     * Test of direct transpose vector solver
     */
    public void testTransVectorSolve() {
        while (true) {
            try {
                Vector b = Matrices.random(A.numRows());
                Vector x = Matrices.random(A.numRows());
                x = A.transSolve(b, x);

                Vector y = A.transMultAdd(-1, x, x.copy().set(b));
                assertEquals(0, y.norm(Vector.Norm.Two), tol);
                assertEquals(Ad, A);
                return;
            } catch (MatrixSingularException e) {
                Utilities.addDiagonal(A, Ad, 1);
            } catch (MatrixNotSPDException e) {
                Utilities.addDiagonal(A, Ad, 1);
            }
        }
    }

    /**
     * Test additions using iterators
     */
    public void testAdd() {
        double alpha = Math.random();
        for (MatrixEntry e : A) {
            A.add(e.row(), e.column(), alpha);
            A.add(e.row(), e.column(), -alpha);
        }
        assertEquals(Ad, A);
    }

    /**
     * Checks that copy is deep, not reference
     */
    public void testCopy() {
        Matrix Ac = A.copy();
        A = A.zero();
        assertEquals(Ad, Ac);
    }

    /**
     * Test iterator get
     */
    public void testIterator() {
        double[][] Ac = new double[A.numRows()][A.numColumns()];
        for (MatrixEntry e : A)
            Ac[e.row()][e.column()] = e.get();
        assertEquals(Ad, Ac);
    }

    /**
     * Test iterator set
     */
    public void testIteratorSet() {
        double alpha = Math.random();
        for (MatrixEntry e : A)
            e.set(e.get() * alpha);
        assertEquals(scale(alpha), A);
    }

    /**
     * Test iterator read and write
     */
    public void testIteratorSetGet() {
        double alpha = Math.random();
        double[][] Ac = new double[A.numRows()][A.numColumns()];
        for (MatrixEntry e : A) {
            Ac[e.row()][e.column()] = e.get();
            e.set(alpha * e.get());
            e.set(e.get() / alpha);
        }
        assertEquals(Ad, Ac);
        assertEquals(Ad, A);
    }

    /**
     * Checks zero()
     */
    public void testZero() {
        assertEquals(zero(), A.zero());
    }

    protected double[][] zero() {
        for (int i = 0; i < Ad.length; ++i)
            for (int j = 0; j < Ad[i].length; ++j)
                Ad[i][j] = 0;
        return Ad;
    }

    /**
     * Cardinality computation
     */
    public void testCardinality() {
        assertEquals(Matrices.cardinality(A), cardinality());
    }

    protected int cardinality() {
        int nz = 0;
        for (int i = 0; i < Ad.length; ++i)
            for (int j = 0; j < Ad[i].length; ++j)
                if (Ad[i][j] != 0.)
                    nz++;
        return nz;
    }

    /**
     * Checks in-place transpose for square matrices
     */
    public void testTransposeInplace() {
        if (A.isSquare())
            assertEquals(transpose(), A.copy().transpose());
    }

    /**
     * Scaling with an arbitrary scalar
     */
    public void testScale() {
        double alpha = Math.random();
        A = A.scale(alpha);
        scale(alpha);
        assertEquals(Ad, A);
    }

    /**
     * Scaling by zero
     */
    public void testZeroScale() {
        A = A.scale(0);
        scale(0);
        assertEquals(Ad, A);
    }

    /**
     * Scaling by one
     */
    public void testOneScale() {
        A = A.scale(1);
        scale(1);
        assertEquals(Ad, A);
    }

    protected double[][] scale(double alpha) {
        for (int i = 0; i < Ad.length; ++i)
            for (int j = 0; j < Ad[i].length; ++j)
                Ad[i][j] *= alpha;
        return Ad;
    }

    /**
     * Checks the 1 norm
     */
    public void testOneNorm() {
        assertEquals(norm1(Ad), A.norm(Matrix.Norm.One), tol);
        assertEquals(Ad, A);
    }

    /**
     * Checks the Frobenius norm
     */
    public void testFrobeniusNorm() {
        assertEquals(normF(Ad), A.norm(Matrix.Norm.Frobenius), tol);
        assertEquals(Ad, A);
    }

    /**
     * Checks the infinity norm
     */
    public void testInfinityNorm() {
        assertEquals(normInf(Ad), A.norm(Matrix.Norm.Infinity), tol);
        assertEquals(Ad, A);
    }

    protected double norm1(double[][] A) {
        double max = 0;
        for (int i = 0; i < A.length; ++i) {
            double rowsum = 0;
            for (int j = 0; j < A[i].length; ++j)
                rowsum += Math.abs(A[i][j]);
            max = Math.max(rowsum, max);
        }
        return max;
    }

    protected double normF(double[][] A) {
        double norm = 0;
        for (int i = 0; i < A.length; ++i)
            for (int j = 0; j < A[i].length; ++j)
                norm += A[i][j] * A[i][j];
        return Math.sqrt(norm);
    }

    protected double normInf(double[][] A) {
        if (A.length == 0)
            return 0;

        double[] columnSum = new double[A[0].length];
        for (int i = 0; i < A.length; ++i)
            for (int j = 0; j < A[i].length; ++j)
                columnSum[j] += Math.abs(A[i][j]);

        double max = 0;
        for (double d : columnSum)
            max = Math.max(max, d);
        return max;
    }

    /**
     * Checks for equality between the matrix and the array
     */
    protected void assertEquals(double[][] Ad, Matrix A) {
        assertTrue(A != null);
        assertTrue(Ad != null);
        assertTrue(A.numRows() == Ad.length);
        for (int i = 0; i < A.numRows(); ++i) {
            assertTrue(A.numColumns() == Ad[i].length);
            for (int j = 0; j < A.numColumns(); ++j)
                assertEquals(Ad[i][j], A.get(i, j), 1e-12);
        }
    }

    /**
     * Checks for equality between two arrays
     */
    protected void assertEquals(double[][] Ad, double[][] Ac) {
        assertTrue(Ac.length == Ad.length);
        for (int i = 0; i < A.numRows(); ++i) {
            assertTrue(Ac[i].length == Ad[i].length);
            for (int j = 0; j < A.numColumns(); ++j)
                assertEquals(Ad[i][j], Ac[i][j], 1e-12);
        }
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

  public static void assertEquals(Matrix expected, Matrix test) {
    assertEquals(expected.numRows(), test.numRows());
    assertEquals(expected.numColumns(), test.numColumns());
    for (int i = 0; i < test.numRows(); i++) {
      for (int j = 0; j < test.numColumns(); j++) {
        assertEquals(expected.get(i, j), test.get(i, j), 0.0001);
      }
    }
  }

  public static void assertEquals(Vector expected, Vector test) {
    assertEquals(expected.size(), test.size());
    for (int i = 0; i < test.size(); i++) {
      assertEquals(expected.get(i), test.get(i), 0.0001);
    }
  }

  public static void assertEqualsOrOpposite(Vector expected, Vector test) {
    assertEquals(expected.size(), test.size());
    for (int i = 0; i < test.size(); i++) {
      double a = expected.get(i);
      double b = test.get(i);
      assertTrue("abs(" + a + ") != abs(" + b + ")", Math.abs(a) - Math.abs(b) < 0.0001);
    }
  }

}
