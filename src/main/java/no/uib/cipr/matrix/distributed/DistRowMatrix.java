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

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;

/**
 * Distributed matrix with row major blocks
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
public class DistRowMatrix extends DistMatrix {

    private static final long serialVersionUID = 3258129167668229176L;

    /**
     * Constructor for DistRowMatrix
     * 
     * @param numRows
     *            Global number of rows
     * @param numColumns
     *            Global number of columns
     * @param comm
     *            Communicator to use
     * @param A
     *            Block diagonal matrix. The sum of the local row sizes of
     *            <code>A</code> must equal the global number, and likewise
     *            with the column sizes.
     * @param B
     *            Off-diagonal matrix part. Its number of columns must equal the
     *            global number of columns, and its number of rows must equal
     *            that of <code>A</code>
     */
    public DistRowMatrix(int numRows, int numColumns, Communicator comm,
            Matrix A, Matrix B) {
        super(numRows, numColumns, comm, A, B);

        if (A.numRows() != B.numRows())
            throw new IllegalArgumentException("A.numRows() != B.numRows()");
        if (B.numColumns() != numColumns)
            throw new IllegalArgumentException("B.numColumns() != numColumns");
    }

    @Override
    public void add(int row, int column, double value) {
        check(row, column);

        if (inA(row, column))
            A.add(row - n[rank], column - m[rank], value);
        else if (local(row, column))
            B.add(row - n[rank], column, value);
        else
            throw new IllegalArgumentException("Row index " + row
                    + " is not local");
    }

    @Override
    public void set(int row, int column, double value) {
        check(row, column);

        if (inA(row, column))
            A.set(row - n[rank], column - m[rank], value);
        else if (local(row, column))
            B.set(row - n[rank], column, value);
        else
            throw new IllegalArgumentException("Row index " + row
                    + " is not local");
    }

    @Override
    public double get(int row, int column) {
        check(row, column);

        if (inA(row, column))
            return A.get(row - n[rank], column - m[rank]);
        else if (local(row, column))
            return B.get(row - n[rank], column);
        else
            throw new IndexOutOfBoundsException("Entry not available locally");
    }

    @Override
    public DistRowMatrix copy() {
        return new DistRowMatrix(numRows, numColumns, comm, A.copy(), B.copy());
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new DistMatrixIterator(n[rank], m[rank], n[rank], 0);
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DistVector && y instanceof DistVector))
            throw new IllegalArgumentException("Vectors must be DistVectors");

        checkMultAdd(x, y);

        DistVector xd = (DistVector) x, yd = (DistVector) y;

        // Recieve the needed components of the global x into the local vector
        scatter.startScatter(xd, locC);

        // Local part
        A.multAdd(alpha, xd.getLocal(), yd.getLocal());

        // Finish communications
        scatter.endSetScatter(xd, locC);

        // Non-local part
        B.multAdd(alpha, locC, yd.getLocal());

        return y;
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DistVector && y instanceof DistVector))
            throw new IllegalArgumentException("Vectors must be DistVectors");

        checkTransMultAdd(x, y);

        // y = 1/alpha * y
        y.scale(1. / alpha);

        // y = A'x + y = A'x + 1/alpha * y

        DistVector xd = (DistVector) x, yd = (DistVector) y;

        // Non-local part
        B.transMult(xd.getLocal(), locR);

        // Send it to the others
        scatter.startGather(locR, yd);

        // Local part
        A.transMultAdd(xd.getLocal(), yd.getLocal());

        // Finish communications, concluding the matrix product
        scatter.endAddGather(locR, yd);

        // y = alpha*y = alpha * A'x + y
        return y.scale(alpha);
    }

    @Override
    public boolean local(int row, int column) {
        return row >= n[rank] && row < n[rank + 1];
    }

    int getRank(int row) {
        int i = 1;
        for (; i < n.length; ++i)
            if (row < n[i])
                break;
        return i - 1;
    }

    @Override
    int[] getDelimiter() {
        return n;
    }

    @Override
    int[] getCommIndices() {
        // Get the unique row indices from B
        Collection<Integer> set = new HashSet<Integer>();
        for (MatrixEntry e : B)
            if (!local(e.column(), e.column()))
                set.add(e.column());

        // Get an array representation
        int[] indices = new int[set.size()];
        int j = 0;
        for (Integer i : set)
            indices[j++] = i;
        return indices;
    }

    @Override
    protected double norm1() {
        // Compute local norm
        double norm = super.norm1();

        // Find the maximum
        double[] recv = new double[1];
        comm.allReduce(new double[] { norm }, recv, Reductions.max());

        return recv[0];
    }

    @Override
    protected double normInf() {
        // Compute as much locally as possible
        double[] columnSum = new double[numColumns];
        for (MatrixEntry e : this)
            columnSum[e.column()] += Math.abs(e.get());

        // Sum in the rest from the other ranks
        double[] recv = new double[numColumns];
        comm.allReduce(columnSum, recv, Reductions.sum());

        // The global maximum
        return max(recv);
    }

    @Override
    public DistRowMatrix zero() {
        super.zero();
        return this;
    }

}
