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

import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * Utilities for the testers
 *
 * NOTE: many of these random matrices are not guaranteed to have solutions
 */
public final class Utilities {

    private Utilities() {
        // No need to instantiate
    }

    /**
     * Populates the sparse matrix in a symmetric fashion
     * 
     * @param A
     *            Matrix to populate
     * @param num
     *            Number of entries on each symmetry band
     * @return The matrix data in dense format
     */
    public static double[][] symmetryPopulate(Matrix A, int num) {
        int n = A.numRows(), m = A.numColumns();
        if (m != n)
            throw new IllegalArgumentException("m != n");
        double[][] values = new double[n][m];
        if (n == 0)
            return values;

        for (int j = 0; j < m; ++j)
            for (int i = 0; i < num; ++i) {
                double value = Math.random();
                int k = (int) (Math.random() * n);
                values[k][j] = value;
                values[j][k] = value;
                A.set(k, j, value);
                A.set(j, k, value);
            }
        return values;
    }

    /**
     * Populates the matrix, putting a given number of entries on each column
     * 
     * @param A
     *            Matrix to populate
     * @param num
     *            Number of entries on each column
     * @return The matrix data in dense format
     */
    public static double[][] columnPopulate(Matrix A, int num) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        if (n == 0)
            return values;

        for (int j = 0; j < m; ++j)
            for (int i = 0; i < num; ++i) {
                double value = Math.random();
                int k = (int) (Math.random() * n);
                values[k][j] = value;
                A.set(k, j, value);
            }
        return values;
    }

    /**
     * Populates the matrix, putting a given number of entries on each row
     * 
     * @param A
     *            Matrix to populate
     * @param num
     *            Number of entries on each row
     * @return The matrix data in dense format
     */
    public static double[][] rowPopulate(Matrix A, int num) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        if (m == 0)
            return values;

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < num; ++j) {
                double value = Math.random();
                int k = (int) (Math.random() * m);
                values[i][k] = value;
                A.set(i, k, value);
            }
        return values;
    }

    /**
     * Gets a row-wise non-zero pattern suitable for creating compressed row
     * matrices
     * 
     * @param n
     *            Number of rows
     * @param m
     *            Number of columns
     * @param b
     *            Number of entries on each row
     */
    public static int[][] getRowPattern(int n, int m, int b) {
        int[][] nz = new int[n][];

        for (int i = 0; i < n; ++i) {
            Set<Integer> row = new HashSet<Integer>();
            for (int j = 0; j < b; ++j)
                row.add(getInt(m));

            nz[i] = new int[row.size()];
            int j = 0;
            for (Integer colind : row)
                nz[i][j++] = colind;
        }

        return nz;
    }

    /**
     * Gets a column-wise non-zero pattern suitable for creating compressed
     * column matrices
     * 
     * @param n
     *            Number of rows
     * @param m
     *            Number of columns
     * @param b
     *            Number of entries on each column
     */
    public static int[][] getColumnPattern(int n, int m, int b) {
        int[][] nz = new int[m][];

        for (int i = 0; i < m; ++i) {
            Set<Integer> column = new HashSet<Integer>();
            for (int j = 0; j < b; ++j)
                column.add(getInt(n));

            nz[i] = new int[column.size()];
            int j = 0;
            for (Integer rowind : column)
                nz[i][j++] = rowind;
        }

        return nz;
    }

    /**
     * Populates the matrix, using the given non-zero pattern
     * 
     * @param A
     *            Matrix to populate
     * @param nz
     *            Column indices on each row
     * @return The matrix data in dense format
     */
    public static double[][] rowPopulate(Matrix A, int[][] nz) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        if (m == 0)
            return values;

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < nz[i].length; ++j) {
                double value = Math.random();
                int k = nz[i][j];
                values[i][k] = value;
                A.set(i, k, value);
            }

        return values;
    }

    /**
     * Populates the matrix, using the given non-zero pattern
     * 
     * @param A
     *            Matrix to populate
     * @param nz
     *            Row indices on each column
     * @return The matrix data in dense format
     */
    public static double[][] columnPopulate(Matrix A, int[][] nz) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        if (n == 0)
            return values;

        for (int j = 0; j < m; ++j)
            for (int i = 0; i < nz[j].length; ++i) {
                double value = Math.random();
                int k = nz[j][i];
                values[k][j] = value;
                A.set(k, j, value);
            }

        return values;
    }

    /**
     * Populates the matrix fully
     * 
     * @param A
     *            Matrix to populate
     * @return The matrix data in dense format
     */
    public static double[][] populate(Matrix A) {
      Random random = new Random();
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) {
              double value = random.nextGaussian();
                values[i][j] = value;
                A.set(i, j, value);
            }
        return values;
    }

    /**
     * Populates the banded matrix
     * 
     * @param A
     *            Matrix to populate
     * @param kl
     *            Number of subdiagonls
     * @param ku
     *            Number of superdiagonals
     * @return The matrix data in dense format
     */
    public static double[][] bandPopulate(Matrix A, int kl, int ku) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                if (j - ku <= i && i <= j + kl) {
                    double value = Math.random();
                    values[i][j] = value;
                    A.set(i, j, value);
                }
        return values;
    }

    /**
     * Populates the banded matrix, but not its main diagonal
     * 
     * @param A
     *            Matrix to populate
     * @param kl
     *            Number of subdiagonls
     * @param ku
     *            Number of superdiagonals
     * @return The matrix data in dense format
     */
    public static double[][] unitBandPopulate(Matrix A, int kl, int ku) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                if (i != j && j - ku <= i && i <= j + kl) {
                    double value = Math.random();
                    values[i][j] = value;
                    A.set(i, j, value);
                }
        return values;
    }

    /**
     * Populates the lower triangular part of the matrix
     * 
     * @param A
     *            Matrix to populate
     * @return The matrix data in dense format
     */
    public static double[][] lowerPopulate(Matrix A) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                if (j <= i) {
                    double value = Math.random();
                    values[i][j] = value;
                    A.set(i, j, value);
                }
        return values;
    }

    /**
     * Populates the upper triangular part of the matrix
     * 
     * @param A
     *            Matrix to populate
     * @return The matrix data in dense format
     */
    public static double[][] upperPopulate(Matrix A) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                if (i <= j) {
                    double value = Math.random();
                    values[i][j] = value;
                    A.set(i, j, value);
                }
        return values;
    }

  public static void upperPopulateGauss(Matrix A) {
    Random random = new Random();
    for (int i = 0; i < A.numRows(); i++)
      for (int j = 0; j <= i; j++)
        A.set(i, j, random.nextGaussian());
  }

    /**
     * Populates the strictly lower triangular part of the matrix
     * 
     * @param A
     *            Matrix to populate
     * @return The matrix data in dense format
     */
    public static double[][] unitLowerPopulate(Matrix A) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                if (j < i) {
                    double value = Math.random();
                    values[i][j] = value;
                    A.set(i, j, value);
                }
        return values;
    }

    /**
     * Populates the strictly upper triangular part of the matrix
     * 
     * @param A
     *            Matrix to populate
     * @return The matrix data in dense format
     */
    public static double[][] unitUpperPopulate(Matrix A) {
        int n = A.numRows(), m = A.numColumns();
        double[][] values = new double[n][m];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                if (i < j) {
                    double value = Math.random();
                    values[i][j] = value;
                    A.set(i, j, value);
                }
        return values;
    }

    /**
     * Puts the upper triangular part into the lower triangle
     */
    public static void lowerSymmetrice(double[][] Ad) {
        int n = Ad.length;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                Ad[j][i] = Ad[i][j];
    }

    /**
     * Puts the lower triangular part into the upper triangle
     */
    public static void upperSymmetrice(double[][] Ad) {
        int n = Ad.length;
        for (int i = 0; i < n; ++i) {
            int m = Ad[0].length;
            for (int j = i + 1; j < m; ++j)
                Ad[j][i] = Ad[i][j];
        }
    }

    /**
     * Sets one on the main diagonal
     */
    public static double[][] unitSet(double[][] Ad) {
        for (int i = 0; i < Ad.length; ++i)
            Ad[i][i] = 1;
        return Ad;
    }

    private final static Random r = new Random();

    /**
     * Returns an integer between zero (inclusive) and max (exclusive)
     */
    public static int getInt(int max) {
        return r.nextInt(max);
    }

    /**
     * Returns an integer between min (inclusive) and max (exclusive)
     */
    public static int getInt(int min, int max) {
        return Math.max(min, getInt(max));
    }

    /**
     * Returns true if the matrix is singular
     */
    public static boolean singular(Matrix A) throws NotConvergedException {
        SVD svd = SVD.factorize(A);
        double[] S = svd.getS();
        for (int i = 0; i < S.length; ++i)
            if (S[i] == 0)
                return true;
        return false;
    }

    /**
     * Returns true if the matrix is positive definite
     */
    public static boolean spd(Matrix A) throws NotConvergedException {
        EVD evd = EVD.factorize(A);
        {
            double[] S = evd.getRealEigenvalues();
            for (int i = 0; i < S.length; ++i)
                if (S[i] <= 0.)
                    return false;
        }
        {
            double[] S = evd.getImaginaryEigenvalues();
            for (int i = 0; i < S.length; ++i)
                if (Math.abs(S[i]) > 1e-10)
                    return false;
        }
        return true;
    }

    /**
     * Populates the given vector, and returns an array containing its values
     */
    public static double[] populate(Vector x) {
        double[] xd = new double[x.size()];
        for (int i = 0; i < x.size(); ++i) {
            double alpha = Math.random();
            xd[i] = alpha;
            x.set(i, alpha);
        }
        return xd;
    }

    /**
     * Populates the given vector, and returns an array containing its values.
     * Only m entries are inserted
     */
    public static double[] populate(Vector x, int m) {
        double[] xd = new double[x.size()];
        for (int i = 0; i < m; ++i) {
            double alpha = Math.random();
            int k = (int) (Math.random() * x.size());
            xd[k] = alpha;
            x.set(k, alpha);
        }
        return xd;
    }

    /**
     * Zeros the given array
     */
    public static void zero(double[][] A) {
        for (int i = 0; i < A.length; ++i)
            Arrays.fill(A[i], 0);
    }

    /**
     * Adds to the diagonal of both the matrix and the array
     */
    public static void addDiagonal(Matrix A, double[][] Ad, double shift) {
        int n = A.numRows(), m = A.numColumns();
        for (int i = 0; i < Math.min(m, n); ++i) {
            A.add(i, i, shift);
            Ad[i][i] += shift;
        }
    }

    /**
     * Adds to the diagonal of the matrix
     */
    public static void addDiagonal(Matrix A, double shift) {
        int n = A.numRows(), m = A.numColumns();
        for (int i = 0; i < Math.min(m, n); ++i)
            A.add(i, i, shift);
    }

}
