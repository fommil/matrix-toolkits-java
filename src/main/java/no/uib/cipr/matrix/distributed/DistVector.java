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

import java.util.Iterator;

import no.uib.cipr.matrix.AbstractVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;

/**
 * Distributed memory vector
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
public class DistVector extends AbstractVector {

    private static final long serialVersionUID = 4048795671438309432L;

    /**
     * Communicator in use
     */
    private Communicator comm;

    /**
     * Local part of the vector
     */
    Vector x;

    /**
     * The subdivisions of the global vector
     */
    int[] n;

    int rank;

    /**
     * Rank and size of the communicator
     */
    private int commSize;

    /**
     * Constructor for DistVector
     * 
     * @param size
     *            Global vector size
     * @param comm
     *            Communicator to use
     * @param x
     *            Local vector, its size cannot exceed the global size, and the
     *            sum of the local vector sizes must equal the global vector
     *            size (this is checked)
     */
    public DistVector(int size, Communicator comm, Vector x) {
        super(size);
        this.comm = comm;
        this.x = x;

        rank = comm.rank();
        commSize = comm.size();
        n = new int[commSize + 1];

        // Find out the sizes of all the parts of the distributed vector
        int[] send = new int[] { x.size() };
        int[][] recv = new int[comm.size()][1];
        comm.allGather(send, recv);

        for (int i = 0; i < commSize; ++i)
            n[i + 1] = n[i] + recv[i][0];

        if (n[commSize] != size)
            throw new IllegalArgumentException("Sum of local vector sizes ("
                    + n[commSize] + ") do not match the global vector size ("
                    + size + ")");
    }

    @Override
    public void set(int index, double value) {
        check(index);

        if (local(index))
            x.set(index - n[rank], value);
        else
            throw new IllegalArgumentException("Index " + index
                    + " is not local");
    }

    @Override
    public void add(int index, double value) {
        check(index);

        if (local(index))
            x.add(index - n[rank], value);
        else
            throw new IllegalArgumentException("Index " + index
                    + " is not local");
    }

    @Override
    public double get(int index) {
        check(index);

        if (local(index))
            return x.get(index - n[rank]);
        else
            throw new IllegalArgumentException("Entry not available locally");
    }

    @Override
    public DistVector copy() {
        return new DistVector(size, comm, x.copy());
    }

    @Override
    public DistVector zero() {
        x.zero();
        return this;
    }

    @Override
    public DistVector scale(double alpha) {
        x.scale(alpha);
        return this;
    }

    @Override
    public DistVector set(double alpha, Vector y) {
        if (!(y instanceof DistVector))
            throw new IllegalArgumentException("Vector must be DistVector");

        checkSize(y);

        Vector yb = ((DistVector) y).getLocal();

        x.set(alpha, yb);

        return this;
    }

    @Override
    public DistVector add(double alpha, Vector y) {
        if (!(y instanceof DistVector))
            throw new IllegalArgumentException("Vector must be DistVector");

        checkSize(y);

        Vector yb = ((DistVector) y).getLocal();

        x.add(alpha, yb);

        return this;
    }

    @Override
    public double dot(Vector y) {
        if (!(y instanceof DistVector))
            throw new IllegalArgumentException("Vector must be a DistVector");

        checkSize(y);

        // Compute local part
        Vector yb = ((DistVector) y).getLocal();
        double ldot = x.dot(yb);

        // Sum all local parts
        double[] recv = new double[1];
        comm.allReduce(new double[] { ldot }, recv, Reductions.sum());

        return recv[0];
    }

    @Override
    protected double norm1() {
        double norm = x.norm(Norm.One);

        double[] recv = new double[1];
        comm.allReduce(new double[] { norm }, recv, Reductions.sum());

        return norm;
    }

    @Override
    protected double norm2_robust() {
        // We'll just call the fast version, as we have to square the norm
        // anyways during communications
        return norm2();
    }

    @Override
    protected double norm2() {
        // Compute local norm squared
        double norm = x.norm(Norm.Two);
        norm *= norm;

        double[] recv = new double[1];
        comm.allReduce(new double[] { norm }, recv, Reductions.sum());

        return Math.sqrt(recv[0]);
    }

    @Override
    protected double normInf() {
        double norm = x.norm(Norm.Infinity);

        double[] recv = new double[1];
        comm.allReduce(new double[] { norm }, recv, Reductions.max());

        return norm;
    }

    /**
     * Returns the local part of the vector
     */
    public Vector getLocal() {
        return x;
    }

    /**
     * Returns which indices are owned by which ranks. The current rank owns the
     * indices <code>n[comm.rank()]</code> (inclusive) to
     * <code>n[comm.rank()+1]</code> (exclusive)
     */
    public int[] getOwnerships() {
        return n;
    }

    /**
     * Returns true if the insertion index is local to this rank, and no
     * communication is needed afterwards.
     */
    public boolean local(int index) {
        return index >= n[rank] && index < n[rank + 1];
    }

    @Override
    public Iterator<VectorEntry> iterator() {
        return new DistVectorIterator();
    }

    /**
     * Gets the communicator associated with this vector
     */
    public Communicator getCommunicator() {
        return comm;
    }

    /**
     * Iterator over a distributed memory vector
     */
    private class DistVectorIterator implements Iterator<VectorEntry> {

        /**
         * Entry of local iterator
         */
        private DistVectorEntry entry;

        /**
         * Iterator of local vector
         */
        private Iterator<VectorEntry> i;

        /**
         * Constructor for DistVectorIterator
         */
        public DistVectorIterator() {
            i = x.iterator();
            entry = new DistVectorEntry();
        }

        public void remove() {
            i.remove();
        }

        public boolean hasNext() {
            return i.hasNext();
        }

        public VectorEntry next() {
            entry.update(i.next());
            return entry;
        }
    }

    /**
     * Entry returned by the iterator
     */
    private class DistVectorEntry implements VectorEntry {

        private int offset;

        private VectorEntry e;

        public DistVectorEntry() {
            offset = n[rank];
        }

        public void update(VectorEntry e) {
            this.e = e;
        }

        public int index() {
            return e.index() + offset;
        }

        public double get() {
            return e.get();
        }

        public void set(double value) {
            e.set(value);
        }
    }
}
