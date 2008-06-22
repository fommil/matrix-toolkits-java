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

package no.uib.cipr.matrix.sparse;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;

/**
 * ILU preconditioner with fill-in. Uses the dual threshold approach of Saad.
 */
public class ILUT implements Preconditioner {

    /**
     * Factorisation matrix
     */
    private final FlexCompRowMatrix LU;

    /**
     * The L and U factors
     */
    private Matrix L, U;

    /**
     * Temporary vector for solving the factorised system
     */
    private final Vector y;

    /**
     * Drop-tolerance
     */
    private final double tau;

    /**
     * Diagonal indices
     */
    private final int[] diagind;

    /**
     * Stores entries in the lower and upper part of the matrix. Used by the
     * dropping rule to determine the largest entries in the two parts of the
     * matrix
     */
    private final List<IntDoubleEntry> lower, upper;

    /**
     * Number of additional entries to keep in the lower and upper part of the
     * factored matrix. The entries of the original matrix are always kept,
     * unless they numerically too small
     */
    private final int p;

    /**
     * Sets up the preconditioner for the given matrix
     * 
     * @param LU
     *            Matrix to use internally. For best performance, its non-zero
     *            pattern should conform to that of the system matrix
     * @param tau
     *            Drop tolerance
     * @param p
     *            Number of entries to keep on each row in of the factored
     *            matrix. This is in addition to the entries of the original
     *            matrix
     */
    public ILUT(FlexCompRowMatrix LU, double tau, int p) {
        if (!LU.isSquare())
            throw new IllegalArgumentException(
                    "ILU only applies to square matrices");

        this.LU = LU;
        this.tau = tau;
        this.p = p;

        int n = LU.numRows();
        lower = new ArrayList<IntDoubleEntry>(n);
        upper = new ArrayList<IntDoubleEntry>(n);
        y = new DenseVector(n);
        diagind = new int[n];
    }

    /**
     * Sets up the preconditioner for the given matrix. Uses a drop-tolerance of
     * 10<sup>-6</sup>, and keeps 50 entries on each row, including the main
     * diagonal and any previous entries in the matrix structure
     * 
     * @param LU
     *            Matrix to use internally. For best performance, its non-zero
     *            pattern should conform to that of the system matrix
     */
    public ILUT(FlexCompRowMatrix LU) {
        this(LU, 1e-6, 25);
    }

    public Vector apply(Vector b, Vector x) {
        // Ly = b, y = L\b
        L.solve(b, y);

        // Ux = L\b = y
        return U.solve(y, x);
    }

    public Vector transApply(Vector b, Vector x) {
        // U'y = b, y = U'\b
        U.transSolve(b, y);

        // L'x = U'\b = y
        return L.transSolve(y, x);
    }

    public void setMatrix(Matrix A) {
        LU.set(A);
        LU.compact();

        factor();
    }

    private void factor() {
        int n = LU.numRows();

        double[] LUi = new double[n];

        // Find the indices to the diagonal entries
        for (int k = 0; k < n; ++k) {
            SparseVector row = LU.getRow(k);
            diagind[k] = findDiagonalIndex(row, k);
            if (diagind[k] < 0)
                throw new RuntimeException("Missing diagonal entry on row "
                        + (k + 1));
        }

        for (int i = 1; i < n; ++i) {

            // Get row i
            SparseVector rowi = LU.getRow(i);

            // Drop tolerance on current row
            double taui = rowi.norm(Vector.Norm.Two) * tau;

            // Store in dense format
            scatter(rowi, LUi);

            for (int k = 0; k < i; ++k) {

                // Get row k
                SparseVector rowk = LU.getRow(k);
                int[] rowIndex = rowk.getIndex();
                int rowUsed = rowk.getUsed();
                double[] rowData = rowk.getData();

                if (rowData[diagind[k]] == 0)
                    throw new RuntimeException("Zero diagonal entry on row "
                            + (k + 1) + " during ILU process");

                double LUik = LUi[k] / rowData[diagind[k]];

                // Check for small elimination entry
                if (Math.abs(LUik) <= taui)
                    continue;

                // Traverse the sparse row k, reducing row i
                for (int j = diagind[k] + 1; j < rowUsed; ++j)
                    LUi[rowIndex[j]] -= LUik * rowData[j];

                // The above has overwritten LUik, so remedy that
                LUi[k] = LUik;
            }

            // Store back into the LU matrix, dropping as needed
            gather(LUi, rowi, taui, i);

            // Update diagonal index on row i if it is outdated
            if (rowi.getIndex()[diagind[i]] != i) {
                diagind[i] = findDiagonalIndex(rowi, i);
                if (diagind[i] < 0)
                    throw new RuntimeException("Missing diagonal entry on row "
                            + (i + 1) + " during ILU process");
            }
        }

        L = new UnitLowerFlexCompRowMatrix(LU, diagind);
        U = new UpperFlexCompRowMatrix(LU, diagind);
    }

    private int findDiagonalIndex(SparseVector v, int k) {
        return no.uib.cipr.matrix.sparse.Arrays.binarySearch(v.getIndex(), k,
                0, v.getUsed());
    }

    /**
     * Copies the sparse vector into a dense array
     */
    private void scatter(SparseVector v, double[] z) {
        int[] index = v.getIndex();
        int used = v.getUsed();
        double[] data = v.getData();
        Arrays.fill(z, 0);
        for (int i = 0; i < used; ++i)
            z[index[i]] = data[i];
    }

    /**
     * Copies the dense array back into the sparse vector, applying a numerical
     * dropping rule and keeping only a given number of entries
     */
    private void gather(double[] z, SparseVector v, double taui, int d) {
        // Number of entries in the lower and upper part of the original matrix
        int nl = 0, nu = 0;
        for (VectorEntry e : v) {
            if (e.index() < d)
                nl++;
            else if (e.index() > d)
                nu++;
        }
        v.zero();

        // Entries in the L part of the vector
        lower.clear();
        for (int i = 0; i < d; ++i)
            if (Math.abs(z[i]) > taui)
                lower.add(new IntDoubleEntry(i, z[i]));

        // Entries in the U part of the vector
        upper.clear();
        for (int i = d + 1; i < z.length; ++i)
            if (Math.abs(z[i]) > taui)
                upper.add(new IntDoubleEntry(i, z[i]));

        // Sort in descending order
        Collections.sort(lower);
        Collections.sort(upper);

        // Always keep the diagonal
        v.set(d, z[d]);

        // Keep at most nl+p lower entries
        for (int i = 0; i < Math.min(nl + p, lower.size()); ++i) {
            IntDoubleEntry e = lower.get(i);
            v.set(e.index, e.value);
        }

        // Keep at most nu+p upper entries
        for (int i = 0; i < Math.min(nu + p, upper.size()); ++i) {
            IntDoubleEntry e = upper.get(i);
            v.set(e.index, e.value);
        }
    }

    /**
     * Stores an integer/value pair, sorted by descending order according to the
     * value
     */
    private static class IntDoubleEntry implements Comparable<IntDoubleEntry> {

        public int index;

        public double value;

        public IntDoubleEntry(int index, double value) {
            this.index = index;
            this.value = value;
        }

        public int compareTo(IntDoubleEntry o) {
            // Descending order, so keep the largest entries first
            if (Math.abs(value) < Math.abs(o.value))
                return 1;
            else if (Math.abs(value) == Math.abs(o.value))
                return 0;
            else
                return -1;
        }

        @Override
        public String toString() {
            return "(" + index + "=" + value + ")";
        }
    }

    /**
     * Unit lower triangular flex-CRS matrix. Only used for triangular solves
     */
    private static class UnitLowerFlexCompRowMatrix extends AbstractMatrix {

        private final FlexCompRowMatrix LU;

        private final int[] diagind;

        public UnitLowerFlexCompRowMatrix(FlexCompRowMatrix LU, int[] diagind) {
            super(LU);
            this.LU = LU;
            this.diagind = diagind;
        }

        @Override
        public Vector solve(Vector b, Vector x) {
            if (!(b instanceof DenseVector) || !(x instanceof DenseVector))
                return super.solve(b, x);

            double[] bd = ((DenseVector) b).getData();
            double[] xd = ((DenseVector) x).getData();

            for (int i = 0; i < numRows; ++i) {

                // Get row i
                SparseVector row = LU.getRow(i);
                int[] index = row.getIndex();
                double[] data = row.getData();

                // xi = bi - sum[j<i] Lij * xj
                double sum = 0;
                for (int j = 0; j < diagind[i]; ++j)
                    sum += data[j] * xd[index[j]];

                xd[i] = bd[i] - sum;
            }

            return x;
        }

        @Override
        public Vector transSolve(Vector b, Vector x) {
            if (!(x instanceof DenseVector))
                return super.transSolve(b, x);

            x.set(b);

            double[] xd = ((DenseVector) x).getData();

            for (int i = numRows - 1; i >= 0; --i) {

                // Get row i
                SparseVector row = LU.getRow(i);
                int[] index = row.getIndex();
                double[] data = row.getData();

                // At this stage, x[i] is known, so move it over to the right
                // hand side for the remaining equations
                for (int j = 0; j < diagind[i]; ++j)
                    xd[index[j]] -= data[j] * xd[i];

            }

            return x;
        }

    }

    /**
     * Upper triangular flex-CRS matrix. Only used for triangular solves
     */
    private static class UpperFlexCompRowMatrix extends AbstractMatrix {

        private final FlexCompRowMatrix LU;

        private final int[] diagind;

        public UpperFlexCompRowMatrix(FlexCompRowMatrix LU, int[] diagind) {
            super(LU);
            this.LU = LU;
            this.diagind = diagind;
        }

        @Override
        public Vector solve(Vector b, Vector x) {
            if (!(b instanceof DenseVector) || !(x instanceof DenseVector))
                return super.solve(b, x);

            double[] bd = ((DenseVector) b).getData();
            double[] xd = ((DenseVector) x).getData();

            for (int i = numRows - 1; i >= 0; --i) {

                // Get row i
                SparseVector row = LU.getRow(i);
                int[] index = row.getIndex();
                int used = row.getUsed();
                double[] data = row.getData();

                // xi = (bi - sum[j>i] Uij * xj) / Uii
                double sum = 0;
                for (int j = diagind[i] + 1; j < used; ++j)
                    sum += data[j] * xd[index[j]];

                xd[i] = (bd[i] - sum) / data[diagind[i]];
            }

            return x;
        }

        @Override
        public Vector transSolve(Vector b, Vector x) {
            if (!(x instanceof DenseVector))
                return super.transSolve(b, x);

            x.set(b);

            double[] xd = ((DenseVector) x).getData();

            for (int i = 0; i < numRows; ++i) {

                // Get row i
                SparseVector row = LU.getRow(i);
                int[] index = row.getIndex();
                int used = row.getUsed();
                double[] data = row.getData();

                // Solve for the current entry
                xd[i] /= data[diagind[i]];

                // Move this known solution over to the right hand side for the
                // remaining equations
                for (int j = diagind[i] + 1; j < used; ++j)
                    xd[index[j]] -= data[j] * xd[i];
            }

            return x;
        }

    }
}
