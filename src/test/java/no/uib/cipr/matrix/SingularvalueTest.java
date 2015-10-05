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

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;
import no.uib.cipr.matrix.JobSVD;
import no.uib.cipr.matrix.sparse.CompDiagMatrix;
import junit.framework.TestCase;

/**
 * Test the singular value solver
 */
public class SingularvalueTest extends TestCase {

    /**
     * Matrix to decompose
     */
    private DenseMatrix A;

    /**
     * Maximum matrix size, to avoid too slow tests
     */
    private final int max = 200;

    public SingularvalueTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
        int n = Utilities.getInt(1, max);
        int m = Utilities.getInt(1, max);
        A = new DenseMatrix(m, n);
        Utilities.populate(A);
    }

    @Override
    protected void tearDown() throws Exception {
        A = null;
    }

    public void testStaticFactorize() throws NotConvergedException {
        assertEquals(A, SVD.factorize(A));
    }

    public void testFactorDefault() throws NotConvergedException {
        SVD svd = new SVD(A.numRows(), A.numColumns());
        assertEquals(A, svd.factor(A.copy()));
    }

    public void testFactorDgesdd() throws NotConvergedException {
        for (JobSVD jobz : JobSVD.values()) {
            SVD svd = new SVD(A.numRows(), A.numColumns(), jobz);
            if (jobz == JobSVD.All) {
                assertEquals(A, svd.factor(A.copy()));
            } else {
                svd.factor(A.copy());
            }
        }
    }

    public void testFactorDgesvd() throws NotConvergedException {

        for (JobSVD jobU : JobSVD.values()) {
            for (JobSVD jobVT : JobSVD.values()) {
                if (jobU != JobSVD.Overwrite || jobVT != JobSVD.Overwrite) {
                    SVD svd = new SVD(A.numRows(), A.numColumns(), jobU, jobVT);

                    if (jobU == JobSVD.All && jobVT == JobSVD.All) {
                        assertEquals(A, svd.factor(A.copy()));
                    } else {
                        svd.factor(A.copy());
                    }
                }
            }
        }
    }

    // Compare between equivalent results of DGESVD and DGESDD
    public void testFactorCompare() throws NotConvergedException {
        for (JobSVD jobz : JobSVD.values()) {
            if (jobz == JobSVD.Overwrite) {
                continue;
            }

            SVD svd1 = new SVD(A.numRows(), A.numColumns(), jobz);
            SVD svd2 = new SVD(A.numRows(), A.numColumns(), jobz, jobz);

            if (jobz == JobSVD.All) {
                assertEquals(svd1.factor(A.copy()), svd2.factor(A.copy()));
            } else if (jobz == JobSVD.Part) {
                assertEquals(svd1.factor(A.copy()), svd2.factor(A.copy()));
            } else if (jobz == JobSVD.None) {
                assertEquals(svd1.factor(A.copy()), svd2.factor(A.copy()));
            }
        }

        SVD svd1 = new SVD(A.numRows(), A.numColumns(), JobSVD.Overwrite);
        SVD svd2;
        if (A.numRows >= A.numColumns) {
            svd2 = new SVD(A.numRows(), A.numColumns(), JobSVD.Overwrite,
                    JobSVD.All);
        } else {
            svd2 = new SVD(A.numRows(), A.numColumns(), JobSVD.All,
                    JobSVD.Overwrite);
        }
        assertEquals(svd1.factor(A.copy()), svd2.factor(A.copy()));
    }

    private void assertEquals(Matrix A, SVD svd) {

        Matrix s = svdToMatrix(svd);

        // Check that A=U*S*Vt
        for (int i = 0; i < A.numRows(); ++i)
            for (int j = 0; j < A.numColumns(); ++j)
                assertEquals(A.get(i, j), s.get(i, j), 1e-12);
    }

    private void assertEquals(SVD svd1, SVD svd2) {
        // Compute U*S*Vt
        Matrix m1 = svdToMatrix(svd1);
        Matrix m2 = svdToMatrix(svd2);

        // Check that A=U*S*Vt
        for (int i = 0; i < m1.numRows(); ++i)
            for (int j = 0; j < m1.numColumns(); ++j)
                assertEquals(m1.get(i, j), m2.get(i, j), 1e-12);
    }

    private DenseMatrix svdToMatrix(SVD svd) {
        if (!svd.hasRightSingularVectors() || !svd.hasLeftSingularVectors()) {
            return new DenseMatrix(0, 0);
        }

        DenseMatrix U = svd.getU();
        DenseMatrix Vt = svd.getVt();
        CompDiagMatrix S = new CompDiagMatrix(U.numColumns, Vt.numRows);

        double[] dia = svd.getS();

        for (int i = 0; i < dia.length && i < S.numColumns && i < S.numRows; i++) {
            S.set(i, i, dia[i]);
        }

        // Compute U*S*Vt
        return (DenseMatrix) U.mult(
                S.mult(Vt, new DenseMatrix(S.numRows(), Vt.numColumns())),
                new DenseMatrix(A.numRows(), A.numColumns()));
    }
}
