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

import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

/**
 * Unit lower triangular CRS matrix. Only used for triangular solves
 */
class UnitLowerCompRowMatrix extends AbstractMatrix {

    private int[] rowptr;

    private int[] colind;

    private double[] data;

    private int[] diagind;

    public UnitLowerCompRowMatrix(CompRowMatrix LU, int[] diagind) {
        super(LU);
        rowptr = LU.getRowPointers();
        colind = LU.getColumnIndices();
        data = LU.getData();
        this.diagind = diagind;
    }

    @Override
    public Vector solve(Vector b, Vector x) {
        if (!(b instanceof DenseVector) || !(x instanceof DenseVector))
            return super.solve(b, x);

        double[] bd = ((DenseVector) b).getData();
        double[] xd = ((DenseVector) x).getData();

        for (int i = 0; i < numRows; ++i) {

            // xi = bi - sum[j<i] Lij * xj
            double sum = 0;
            for (int j = rowptr[i]; j < diagind[i]; ++j)
                sum += data[j] * xd[colind[j]];

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

        for (int i = numRows - 1; i >= 0; --i)

            // At this stage, x[i] is known, so move it over to the right hand
            // side for the remaining equations
            for (int j = rowptr[i]; j < diagind[i]; ++j)
                xd[colind[j]] -= data[j] * xd[i];

        return x;
    }

}
