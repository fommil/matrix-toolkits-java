/*
 * Copyright (C) 2015 Rogério Pontes
 *
 */

package no.uib.cipr.matrix;

import org.junit.Test;

/**
 * Test of the Khatri Rao Multiplication.
 */
public class KhatriRaoTest {

    @Test
    public void testEqualSizeKhatriRao() {
        Matrix A = new DenseMatrix(new double[][]{{1, 2, 3}, {4, 5, 6},
                {7, 8, 9}});
        Matrix B = new DenseMatrix(new double[][]{{1, 4, 7}, {2, 5, 8},
                {3, 6, 9}});
        Matrix res = new DenseMatrix(new double[][]{{1, 8, 21}, {2, 10, 24},
                {3, 12, 27}, {4, 20, 42}, {8, 25, 48}, {12, 30, 54},
                {7, 32, 63}, {14, 40, 72}, {21, 48, 81}});
        KR mult = new KR(A, B);
        Matrix C = new DenseMatrix(A.numRows() * B.numRows(), A.numColumns());
        C = mult.multiply(C);
        MatrixTestAbstract.assertMatrixEquals(C, res);

    }

    @Test
    public void testEqualColumnKhatriRao() {
        Matrix A = new DenseMatrix(new double[][]{{1, 2}, {4, 2}, {7, 8}});
        Matrix B = new DenseMatrix(new double[][]{{1, 4}, {8, 5}, {5, 4},
                {10, 20}, {30, 40}});
        Matrix res = new DenseMatrix(new double[][]{{1, 8}, {8, 10}, {5, 8},
                {10, 40}, {30, 80}, {4, 8}, {32, 10}, {20, 8}, {40, 40},
                {120, 80}, {7, 32}, {56, 40}, {35, 32}, {70, 160}, {210, 320}});
        KR mult = new KR(A, B);
        Matrix C = new DenseMatrix(A.numRows() * B.numRows(), A.numColumns());
        C = mult.multiply(C);
        MatrixTestAbstract.assertMatrixEquals(C, res);

    }

    @Test
    public void testEqualRowKhatriRao() {
        Matrix A = new DenseMatrix(new double[][]{{1, 2, 4, 5, 6},
                {2, 8, 9, 2, 3}});
        Matrix B = new DenseMatrix(new double[][]{{10, 20, 6, 12, 10},
                {20, 5, 6, 7, 10}});

        Matrix res = new DenseMatrix(new double[][]{{10, 40, 24, 60, 60},
                {20, 10, 24, 35, 60}, {20, 160, 54, 24, 30},
                {40, 40, 54, 14, 30}});
        KR mult = new KR(A, B);
        Matrix C = new DenseMatrix(A.numRows() * B.numRows(), A.numColumns());
        C = mult.multiply(C);
        MatrixTestAbstract.assertMatrixEquals(C, res);
    }

}
