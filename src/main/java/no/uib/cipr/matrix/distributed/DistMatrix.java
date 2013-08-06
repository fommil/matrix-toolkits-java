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
import java.util.Iterator;

import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.distributed.SuperIterator.SuperIteratorEntry;

/**
 * Distributed memory matrix
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
abstract class DistMatrix extends AbstractMatrix {

    /**
     * Communicator in use
     */
    final Communicator comm;

    /**
     * Block diagonal part
     */
    final Matrix A;

    /**
     * Off-diagonal part
     */
    final Matrix B;

    /**
     * Offsets into the local matrix
     */
    final int[] n, m;

    /**
     * Rank and size of the communicator
     */
    final int rank, size;

    /**
     * Local vector caches, for scatter/gather operations. The first with size
     * equal numRows, the other of numColumns size
     */
    final Vector locR, locC;

    /**
     * Scatters global vectors into local, and the other way around
     */
    final VectorScatter scatter;

    /**
     * Constructor for DistMatrix
     */
    public DistMatrix(int numRows, int numColumns, Communicator comm, Matrix A,
            Matrix B) {
        super(numRows, numColumns);
        this.comm = comm;
        this.A = A;
        this.B = B;

        locR = new DenseVector(numRows);
        locC = new DenseVector(numColumns);
        rank = comm.rank();
        size = comm.size();
        n = new int[size + 1];
        m = new int[size + 1];

        // Find out the sizes of all the parts of the distributed vector
        int[] send = new int[] { A.numRows(), A.numColumns() };
        int[][] recv = new int[size][2];
        comm.allGather(send, recv);

        for (int i = 0; i < size; ++i) {
            n[i + 1] = n[i] + recv[i][0]; // rows
            m[i + 1] = m[i] + recv[i][1]; // columns
        }

        if (n[size] != numRows)
            throw new IllegalArgumentException("Sum of local row sizes ("
                    + n[size] + ") do not match the global row size ("
                    + numRows + ")");
        if (m[size] != numColumns)
            throw new IllegalArgumentException("Sum of local column sizes ("
                    + m[size] + ") do not match the global column size ("
                    + numColumns + ")");

        scatter = scatterSetup();
    }

    /**
     * Sets up vector scatter for matrix/vector products. Collective operation
     */
    private VectorScatter scatterSetup() {
        // Get all the indices needing communication
        int[] ind = getCommIndices();
        Arrays.sort(ind);

        // Find who owns what
        int[] N = getDelimiter();

        // Then get the number of entries to recieve and their indices
        int[][] recv = new int[size][1];
        int[][] recvI = new int[size][];
        for (int k = 0, i = 0; k < size; ++k) {

            // Number of entries
            for (int l = 0, I = i; I < ind.length && ind[I] < N[k + 1]; ++I, ++l)
                recv[k][0]++;

            // The indices
            recvI[k] = new int[recv[k][0]];
            for (int l = 0, I = i; I < ind.length && ind[I] < N[k + 1]; ++I, ++l)
                recvI[k][l] = ind[I];

            i += recv[k][0];
        }

        // Get the number of entries to send
        int[][] send = new int[size][1];
        comm.allToAll(recv, send);

        // Then the indices to send
        int[][] sendI = new int[size][];
        for (int i = 0; i < size; ++i)
            sendI[i] = new int[send[i][0]];

        // Send what we need to recieve, and in return, get what the other
        // threads need
        comm.allToAll(recvI, sendI);

        // Create vector scatter object
        return new VectorScatter(comm, sendI, recvI);
    }

    /**
     * Returns delimiters
     */
    abstract int[] getDelimiter();

    /**
     * Returns indices needing communication (off the block diagonal)
     */
    abstract int[] getCommIndices();

    /**
     * Returns which rows are owned by which ranks. The current rank owns the
     * rows <code>n[comm.rank()]</code> (inclusive) to
     * <code>n[comm.rank()+1]</code> (exclusive)
     */
    public int[] getRowOwnerships() {
        return n;
    }

    /**
     * Returns which columns are owned by which ranks. The current rank owns the
     * columns <code>m[comm.rank()]</code> (inclusive) to
     * <code>m[comm.rank()+1]</code> (exclusive)
     */
    public int[] getColumnOwnerships() {
        return m;
    }

    /**
     * Returns the diagonal block matrix
     */
    public Matrix getBlock() {
        return A;
    }

    /**
     * Returns the off-diagonal matrix
     */
    public Matrix getOff() {
        return B;
    }

    @Override
    public DistMatrix zero() {
        A.zero();
        B.zero();
        return this;
    }

    @Override
    protected double max() {
        // Compute local norms
        double normA = A.norm(Norm.Maxvalue), normB = B.norm(Norm.Maxvalue);

        // Find global maximum
        double[] recv = new double[2];
        comm.allReduce(new double[] { normA, normB }, recv, Reductions.max());

        return recv[0] + recv[1];
    }

    @Override
    protected double normF() {
        // Compute local norms
        double normA = A.norm(Norm.Frobenius), normB = B.norm(Norm.Frobenius);
        normA *= normA;
        normB *= normB;

        // Sum the global norms
        double[] recv = new double[2];
        comm.allReduce(new double[] { normA, normB }, recv, Reductions.sum());

        return Math.sqrt(recv[0] + recv[1]);
    }

    /**
     * Returns true if the insertion indices are local to this rank, and no
     * communication is required afterwards. However, you still need to call
     * <code>flushAssembly</code> to set up things like matrix/vector
     * multiplication
     */
    public abstract boolean local(int row, int column);

    boolean inA(int row, int column) {
        return row >= n[rank] && row < n[rank + 1] && column >= m[rank]
                && column < m[rank + 1];
    }

    @Override
    public Matrix rank1(double alpha, Vector x, Vector y) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Matrix rank2(double alpha, Vector x, Vector y) {
        throw new UnsupportedOperationException();
    }

    /**
     * Gets the communicator associated with this matrix
     */
    public Communicator getCommunicator() {
        return comm;
    }

    /**
     * Iterator for a distributed memory matrix
     */
    class DistMatrixIterator implements Iterator<MatrixEntry> {

        /**
         * Iterates over each column vector
         */
        private SuperIterator<Matrix, MatrixEntry> iterator;

        /**
         * Entry returned
         */
        private DistMatrixEntry entry;

        private int rowAOffset, columnAOffset, rowBOffset, columnBOffset;

        public DistMatrixIterator(int rowAOffset, int columnAOffset,
                int rowBOffset, int columnBOffset) {
            this.rowAOffset = rowAOffset;
            this.rowBOffset = rowBOffset;
            this.columnAOffset = columnAOffset;
            this.columnBOffset = columnBOffset;

            iterator = new SuperIterator<Matrix, MatrixEntry>(Arrays.asList(A,
                    B));
            entry = new DistMatrixEntry();
        }

        public boolean hasNext() {
            return iterator.hasNext();
        }

        public MatrixEntry next() {
            SuperIteratorEntry<MatrixEntry> se = iterator.next();
            if (se.index() == 0)
                // Block diagonal part
                entry.update(rowAOffset, columnAOffset, se.get());
            else
                // Off diagonal block
                entry.update(rowBOffset, columnBOffset, se.get());
            return entry;
        }

        public void remove() {
            iterator.remove();
        }
    }

    /**
     * Entry of a distributed memory matrix
     */
    private static class DistMatrixEntry implements MatrixEntry {

        private int row, column;

        private MatrixEntry entry;

        public void update(int rowOffset, int columnOffset, MatrixEntry entry) {
            row = rowOffset + entry.row();
            column = columnOffset + entry.column();
            this.entry = entry;
        }

        public int row() {
            return row;
        }

        public int column() {
            return column;
        }

        public double get() {
            return entry.get();
        }

        public void set(double value) {
            entry.set(value);
        }

    }

}
