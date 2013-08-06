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

package no.uib.cipr.matrix.distributed;

import java.util.Arrays;

import no.uib.cipr.matrix.DenseLU;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;
import no.uib.cipr.matrix.sparse.Preconditioner;

/**
 * Two level preconditioner. Uses a block preconditioner as a subdomain solver,
 * and algebraically constructs a coarse grid correcion operator
 *
 * @deprecated the <code>no.uib.cipr.matrix.distributed</code> package has been deprecated because
 * of a number of hard to fix concurrency bugs. It is distributed only for backwards compatibility,
 * but is not recommended. The utility of this package is questionable, as it does not allow
 * distribution of computation between JVMs or across a network. For many people, distributed
 * computing of multiple matrices can be achieved at a user-level through the
 * <a href="http://jppf.org">JPPF Framework</a>.
 * Users who need to deal with few very large matrices may wish to implement their own storage classes
 * and solvers using JPPF, but this will not be supported directly in matrix-toolkits-java.
 */
@Deprecated
public class TwoLevelPreconditioner extends BlockDiagonalPreconditioner {

    private final static int root = 0;

    private final DistMatrix A;

    private final Communicator comm;

    private final int rank, size;

    private final int[] indexToRank;

    private final DistVector z, r;

    private final DenseMatrix A0;

    private final DenseVector b0;

    private final DenseLU lu;

    private final double[] Ai;

    private final double[][] Ai0;

    private final boolean row;

    private final double[][] zi0;

    public TwoLevelPreconditioner(Preconditioner prec, DistRowMatrix A,
            DistVector z) {
        super(prec);
        this.A = A;
        this.z = z;
        this.r = z.copy();
        row = true;

        indexToRank = createIndexToRank(A.numColumns(), A.getColumnOwnerships());

        this.comm = A.getCommunicator();
        rank = comm.rank();
        size = comm.size();

        A0 = new DenseMatrix(size, size);
        b0 = new DenseVector(size);
        lu = new DenseLU(size, size);

        Ai = new double[size];
        if (rank == root) {
            Ai0 = new double[size][size];
            zi0 = new double[size][1];
        } else {
            Ai0 = null;
            zi0 = null;
        }
    }

    public TwoLevelPreconditioner(Preconditioner prec, DistColMatrix A,
            DistVector z) {
        super(prec);
        this.A = A;
        this.z = z;
        this.r = z.copy();
        row = false;

        indexToRank = createIndexToRank(A.numColumns(), A.getColumnOwnerships());

        this.comm = A.getCommunicator();
        rank = comm.rank();
        size = comm.size();

        A0 = new DenseMatrix(size, size);
        b0 = new DenseVector(size);
        lu = new DenseLU(size, size);

        Ai = new double[size];
        if (rank == root) {
            Ai0 = new double[size][size];
            zi0 = new double[size][1];
        } else {
            Ai0 = null;
            zi0 = null;
        }
    }

    /**
     * Creates a mapping from a row/column index to the associated rank
     * 
     * @param size
     *            Number of rows/columns
     * @param n
     *            Row/column ownerships
     */
    private int[] createIndexToRank(int size, int[] n) {
        int[] indexToRank = new int[size];

        for (int i = 0; i < n.length - 1; ++i)
            for (int j = n[i]; j < n[i + 1]; ++j)
                indexToRank[j] = i;

        return indexToRank;
    }

    @Override
    public Vector apply(Vector b, Vector x) {
        if (!(b instanceof DistVector) || !(x instanceof DistVector))
            throw new IllegalArgumentException("Vectors must be DistVectors");

        boolean transpose = false;

        return apply(b, x, transpose);
    }

    @Override
    public Vector transApply(Vector b, Vector x) {
        if (!(b instanceof DistVector) || !(x instanceof DistVector))
            throw new IllegalArgumentException("Vectors must be DistVectors");

        boolean transpose = true;

        return apply(b, x, transpose);
    }

    private Vector apply(Vector b, Vector x, boolean transpose) {
        // Calculate R * (b - A*x)
        calculateCoarseResidual(b, x, transpose);

        // Calculate A0 \ R * (b - A*x)
        if (rank == root)
            solveCoarseSystem(transpose);

        // Update x += R^T A0 \ R * (b - A*x)
        updateWithCoarseCorrection(x);

        return applyBlockPreconditioner(b, x, transpose);
    }

    /**
     * Calculates the residual, and projects it onto the coarse grid
     */
    private void calculateCoarseResidual(Vector b, Vector x, boolean transpose) {
        // z = b - A*x
        if (transpose)
            A.transMultAdd(-1, x, z.set(b));
        else
            A.multAdd(-1, x, z.set(b));

        // Piecewise constant projection R * (b - A*x)
        double zi = 0;
        for (VectorEntry e : z.getLocal())
            zi += e.get();

        // Gather the results on rank zero
        double[] zij = new double[] { zi };
        comm.gather(zij, zi0, root);
    }

    /**
     * Solves the coarse system. Only the root thread should call
     */
    private void solveCoarseSystem(boolean transpose) {
        for (int i = 0; i < size; ++i)
            b0.set(i, zi0[i][0]);

        if (transpose)
            lu.transSolve(new DenseMatrix(b0, false));
        else
            lu.solve(new DenseMatrix(b0, false));

        double[] data = b0.getData();
        for (int i = 0; i < size; ++i)
            zi0[i][0] = data[i];
    }

    /**
     * Updates the global solution with the coarse correction
     */
    private void updateWithCoarseCorrection(Vector x) {
        DistVector xd = (DistVector) x;

        double[] zij = new double[1];
        comm.scatter(zi0, zij, root);
        for (VectorEntry e : xd.getLocal())
            e.set(e.get() + zij[0]);
    }

    /**
     * Applies the local block preconditioner
     */
    private Vector applyBlockPreconditioner(Vector b, Vector x,
            boolean transpose) {
        // Recalculate residual
        if (transpose)
            A.transMultAdd(-1, x, z.set(b));
        else
            A.multAdd(-1, x, z.set(b));

        // Apply block preconditioner on the residual
        r.set(b);
        if (transpose)
            super.transApply(z, r);
        else
            super.apply(z, r);

        return x.add(r);
    }

    @Override
    public void setMatrix(Matrix A) {
        if (!(A instanceof DistMatrix))
            throw new IllegalArgumentException(
                    "A is not a DistRowMatrix or a DistColMatrix");

        Matrix Ad = this.A.getBlock();
        Matrix Ao = this.A.getOff();

        super.setMatrix(A);

        Arrays.fill(Ai, 0);
        for (MatrixEntry e : Ad)
            Ai[rank] += e.get();

        if (row) {
            for (MatrixEntry e : Ao)
                Ai[indexToRank[e.column()]] += e.get();

            comm.gather(Ai, Ai0, root);
            if (rank == root)
                for (int i = 0; i < size; ++i)
                    for (int j = 0; j < size; ++j)
                        A0.set(i, j, Ai0[i][j]);
        } else {
            for (MatrixEntry e : Ao)
                Ai[indexToRank[e.row()]] += e.get();

            comm.gather(Ai, Ai0, root);
            if (rank == root)
                for (int i = 0; i < size; ++i)
                    for (int j = 0; j < size; ++j)
                        A0.set(i, j, Ai0[j][i]);
        }

        if (rank == root)
            lu.factor(A0);
    }
}
