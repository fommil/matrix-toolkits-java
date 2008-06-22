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

import java.util.Arrays;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 * Incomplete Cholesky preconditioner without fill-in using a compressed row
 * matrix as internal storage
 */
public class ICC implements Preconditioner {

    /**
     * Factorisation matrix
     */
    private final CompRowMatrix R;

    /**
     * Triangular view onto R for solution purposes
     */
    private Matrix Rt;

    /**
     * Temporary vector for solving the factorised system
     */
    private final Vector y;

    /**
     * Sets up the ICC preconditioner
     * 
     * @param R
     *            Matrix to use internally. For best performance, its non-zero
     *            pattern must conform to that of the system matrix
     */
    public ICC(CompRowMatrix R) {
        if (!R.isSquare())
            throw new IllegalArgumentException(
                    "ICC only applies to square matrices");

        this.R = R;
        int n = R.numRows();
        y = new DenseVector(n);
    }

    public Vector apply(Vector b, Vector x) {
        // R'y = b, y = R'\b
        Rt.transSolve(b, y);

        // Rx = R'\b = y
        return Rt.solve(y, x);
    }

    public Vector transApply(Vector b, Vector x) {
        return apply(b, x);
    }

    public void setMatrix(Matrix A) {
        R.set(A);

        factor();
    }

    private void factor() {
        int n = R.numRows();

        // Internal CRS matrix storage
        int[] colind = R.getColumnIndices();
        int[] rowptr = R.getRowPointers();
        double[] data = R.getData();

        // Temporary storage of a dense row
        double[] Rk = new double[n];

        // Find the indices to the diagonal entries
        int[] diagind = findDiagonalIndices(n, colind, rowptr);

        // Go down along the main diagonal
        for (int k = 0; k < n; ++k) {

            // Expand current row to dense storage
            Arrays.fill(Rk, 0);
            for (int i = rowptr[k]; i < rowptr[k + 1]; ++i)
                Rk[colind[i]] = data[i];

            for (int i = 0; i < k; ++i) {

                // Get the current diagonal entry
                double Rii = data[diagind[i]];

                if (Rii == 0)
                    throw new RuntimeException("Zero pivot encountered on row "
                            + (i + 1) + " during ICC process");

                // Elimination factor
                double Rki = Rk[i] / Rii;

                if (Rki == 0)
                    continue;

                // Traverse the sparse row i, reducing on row k
                for (int j = diagind[i] + 1; j < rowptr[i + 1]; ++j)
                    Rk[colind[j]] -= Rki * data[j];
            }

            // Store the row back into the factorisation matrix
            if (Rk[k] == 0)
                throw new RuntimeException(
                        "Zero diagonal entry encountered on row " + (k + 1)
                                + " during ICC process");
            double sqRkk = Math.sqrt(Rk[k]);

            for (int i = diagind[k]; i < rowptr[k + 1]; ++i)
                data[i] = Rk[colind[i]] / sqRkk;
        }

        Rt = new UpperCompRowMatrix(R, diagind);
    }

    private int[] findDiagonalIndices(int m, int[] colind, int[] rowptr) {
        int[] diagind = new int[m];

        for (int k = 0; k < m; ++k) {
            diagind[k] = no.uib.cipr.matrix.sparse.Arrays.binarySearch(colind,
                    k, rowptr[k], rowptr[k + 1]);

            if (diagind[k] < 0)
                throw new RuntimeException("Missing diagonal entry on row "
                        + (k + 1));
        }

        return diagind;
    }
}
