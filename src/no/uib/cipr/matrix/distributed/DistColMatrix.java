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
 * Distributed matrix with column major blocks
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
public class DistColMatrix extends DistMatrix {

    private static final long serialVersionUID = 3618421523053359673L;

    /**
     * Constructor for DistColMatrix
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
     *            Off-diagonal matrix part. Its number of rows must equal the
     *            global number of rows, and its number of columns must equal
     *            that of <code>A</code>
     */
    public DistColMatrix(int numRows, int numColumns, Communicator comm,
            Matrix A, Matrix B) {
        super(numRows, numColumns, comm, A, B);

        if (A.numColumns() != B.numColumns())
            throw new IllegalArgumentException(
                    "A.numColumns() != B.numColumns()");
        if (B.numRows() != numRows)
            throw new IllegalArgumentException("B.numRows() != numRows");
    }

    @Override
    public void add(int row, int column, double value) {
        check(row, column);

        if (inA(row, column))
            A.add(row - n[rank], column - m[rank], value);
        else if (local(row, column))
            B.add(row, column - m[rank], value);
        else
            throw new IllegalArgumentException("Column index " + column
                    + " is not local");
    }

    @Override
    public void set(int row, int column, double value) {
        check(row, column);

        if (inA(row, column))
            A.set(row - n[rank], column - m[rank], value);
        else if (local(row, column))
            B.set(row, column - m[rank], value);
        else
            throw new IllegalArgumentException("Column index " + column
                    + " is not local");
    }

    @Override
    public double get(int row, int column) {
        check(row, column);

        if (inA(row, column))
            return A.get(row - n[rank], column - m[rank]);
        else if (local(row, column))
            return B.get(row, column - m[rank]);
        else
            throw new IndexOutOfBoundsException("Entry not available locally");
    }

    @Override
    public DistColMatrix copy() {
        return new DistColMatrix(numRows, numColumns, comm, A.copy(), B.copy());
    }

    @Override
    public Iterator<MatrixEntry> iterator() {
        return new DistMatrixIterator(n[rank], m[rank], 0, m[rank]);
    }

    @Override
    public Vector multAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DistVector && y instanceof DistVector))
            throw new IllegalArgumentException("Vectors must be DistVectors");

        checkMultAdd(x, y);

        // y = 1/alpha * y
        y.scale(1. / alpha);

        // y = A*x + y = A*x + 1/alpha * y

        DistVector xd = (DistVector) x, yd = (DistVector) y;

        // Non-local part
        B.mult(xd.getLocal(), locR);

        // Send it to the others
        scatter.startGather(locR, yd);

        // Local part
        A.multAdd(xd.getLocal(), yd.getLocal());

        // Finish communications, concluding the matrix product
        scatter.endAddGather(locR, yd);

        // y = alpha*y = alpha * A*x + y
        return y.scale(alpha);
    }

    @Override
    public Vector transMultAdd(double alpha, Vector x, Vector y) {
        if (!(x instanceof DistVector && y instanceof DistVector))
            throw new IllegalArgumentException("Vectors must be DistVectors");

        checkTransMultAdd(x, y);

        DistVector xd = (DistVector) x, yd = (DistVector) y;

        // Recieve the needed components of the global x into the local vector
        scatter.startScatter(xd, locC);

        // Local part
        A.transMultAdd(alpha, xd.getLocal(), yd.getLocal());

        // Finish communications
        scatter.endSetScatter(xd, locC);

        // Non-local part
        B.transMultAdd(alpha, locC, yd.getLocal());

        return y;
    }

    @Override
    public boolean local(int row, int column) {
        return column >= m[rank] && column < m[rank + 1];
    }

    int getRank(int column) {
        int i = 1;
        for (; i < m.length; ++i)
            if (column < m[i])
                break;
        return i - 1;
    }

    @Override
    int[] getDelimiter() {
        return m;
    }

    @Override
    int[] getCommIndices() {
        // Get the unique row indices from B
        Collection<Integer> set = new HashSet<Integer>();
        for (MatrixEntry e : B)
            if (!local(e.row(), e.row()))
                set.add(e.row());

        // Get an array representation
        int[] indices = new int[set.size()];
        int j = 0;
        for (Integer i : set)
            indices[j++] = i;
        return indices;
    }

    @Override
    protected double norm1() {
        // Compute as much locally as possible
        double[] rowSum = new double[numRows];
        for (MatrixEntry e : this)
            rowSum[e.row()] += Math.abs(e.get());

        // Sum in the rest from the other ranks
        double[] recv = new double[numRows];
        comm.allReduce(rowSum, recv, Reductions.sum());

        // The global maximum
        return max(recv);
    }

    @Override
    protected double normInf() {
        // Compute local norm
        double norm = super.normInf();

        // Find the maximum
        double[] recv = new double[1];
        comm.allReduce(new double[] { norm }, recv, Reductions.max());

        return recv[0];
    }

    @Override
    public DistColMatrix zero() {
        super.zero();
        return this;
    }

}
