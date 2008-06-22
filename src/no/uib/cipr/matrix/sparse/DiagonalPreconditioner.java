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

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 * Diagonal preconditioner. Uses the inverse of the diagonal as preconditioner
 */
public class DiagonalPreconditioner implements Preconditioner {

    /**
     * This contains the inverse of the diagonal
     */
    private double[] invdiag;

    /**
     * Constructor for DiagonalPreconditioner
     * 
     * @param n
     *            Problem size (number of rows)
     */
    public DiagonalPreconditioner(int n) {
        invdiag = new double[n];
    }

    public Vector apply(Vector b, Vector x) {
        if (!(x instanceof DenseVector) || !(b instanceof DenseVector))
            throw new IllegalArgumentException("Vector must be DenseVectors");

        double[] xd = ((DenseVector) x).getData();
        double[] bd = ((DenseVector) b).getData();

        for (int i = 0; i < invdiag.length; ++i)
            xd[i] = bd[i] * invdiag[i];

        return x;
    }

    public Vector transApply(Vector b, Vector x) {
        return apply(b, x);
    }

    public void setMatrix(Matrix A) {
        if (A.numRows() != invdiag.length)
            throw new IllegalArgumentException(
                    "Matrix size differs from preconditioner size");

        for (int i = 0; i < invdiag.length; ++i) {
            invdiag[i] = A.get(i, i);
            if (invdiag[i] == 0) // Avoid zero-division
                throw new RuntimeException("Zero diagonal on row " + (i + 1));
            else
                invdiag[i] = 1 / invdiag[i];
        }
    }

}
