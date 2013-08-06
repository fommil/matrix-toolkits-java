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

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.Exchanger;

import no.uib.cipr.matrix.distributed.Communicator.SendRecv;

/**
 * Collective communications between all the participating threads.
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
public class CollectiveCommunications {

    /**
     * Size of the whole collective
     */
    final int size;

    /**
     * Exchanges pipes between the threads. Only used for setup of the
     * point-to-point communications
     */
    private final List<List<Exchanger<SendRecv>>> ex;

    /**
     * Sets up a collective of the given size
     * 
     * @param size
     *            Number of members of the collective
     */
    public CollectiveCommunications(int size) {
        if (size < 1)
            throw new IllegalArgumentException("size < 1");

        this.size = size;

        barrier = new CyclicBarrier(size);
        broadcast = new Broadcast();
        gather = new Gather();
        scatter = new Scatter();
        allGather = new AllGather();
        allToAll = new AllToAll();
        reduce = new Reduce();

        ex = new ArrayList<List<Exchanger<SendRecv>>>();
        for (int i = 0; i < size; ++i) {
            List<Exchanger<SendRecv>> iex = new ArrayList<Exchanger<SendRecv>>();

            // Get from symmetry
            for (int j = 0; j < i; ++j)
                iex.add(ex.get(j).get(i));

            // Main diagonal (ignored)
            iex.add(null);

            // New entries
            for (int j = i + 1; j < size; ++j)
                iex.add(new Exchanger<SendRecv>());

            ex.add(iex);
        }
    }

    /**
     * Gets the size of the collective
     */
    public int size() {
        return size;
    }

    /**
     * Creates a communicator for point-to-point data-exchange. This method
     * should be called by all the threads in the collective
     * 
     * @param rank
     *            Rank of the communicator
     */
    public Communicator createCommunicator(int rank) {
        if (rank < 0 || rank >= size)
            throw new IllegalArgumentException("rank < 0 || rank >= size");

        return new Communicator(rank, ex.get(rank), this);
    }

    static void await(CyclicBarrier barrier) {
        try {
            barrier.await();
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        } catch (BrokenBarrierException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Blocks the caller until all threads have called
     */
    void barrier() {
        await(barrier);
    }

    private final CyclicBarrier barrier;

    void broadcast(Object buffer, int root, int rank) {
        broadcast.buffer[rank] = buffer;
        if (rank == root)
            broadcast.root = root;

        await(broadcast.barrier);
    }

    private final Broadcast broadcast;

    private class Broadcast implements Runnable {

        final CyclicBarrier barrier = new CyclicBarrier(size, this);

        int root;

        final Object[] buffer = new Object[size];

        public void run() {
            int length = Array.getLength(buffer[root]);
            for (int i = 0; i < size; ++i)
                System.arraycopy(buffer[root], 0, buffer[i], 0, length);
        }
    }

    void gather(Object sendbuf, Object[] recvbuf, int root, int rank) {
        gather.setSend(sendbuf, rank);
        if (rank == root)
            gather.recvbuf = recvbuf;

        await(gather.barrier);
    }

    private final Gather gather;

    private class Gather implements Runnable {

        CyclicBarrier barrier = new CyclicBarrier(size, this);

        Object[] recvbuf = new Object[size], sendbuf = new Object[size];

        private final int[] length = new int[size];

        public void setSend(Object sendbuf, int rank) {
            this.sendbuf[rank] = sendbuf;
            length[rank] = Array.getLength(sendbuf);
        }

        public void run() {
            for (int i = 0; i < size; ++i)
                System.arraycopy(sendbuf[i], 0, recvbuf[i], 0, length[i]);
        }
    }

    void scatter(Object[] sendbuf, Object recvbuf, int root, int rank) {
        scatter.setRecv(recvbuf, rank);
        if (rank == root)
            scatter.sendbuf = sendbuf;

        await(scatter.barrier);
    }

    private final Scatter scatter;

    private class Scatter implements Runnable {

        final CyclicBarrier barrier = new CyclicBarrier(size, this);

        Object[] sendbuf = new Object[size], recvbuf = new Object[size];

        private final int[] length = new int[size];

        public void setRecv(Object recvbuf, int rank) {
            this.recvbuf[rank] = recvbuf;
            length[rank] = Array.getLength(recvbuf);
        }

        public void run() {
            for (int i = 0; i < size; ++i)
                System.arraycopy(sendbuf[i], 0, recvbuf[i], 0, length[i]);
        }
    }

    void allGather(Object sendbuf, Object[] recvbuf, int rank) {
        allGather.setSendRecv(sendbuf, recvbuf, rank);

        await(allGather.barrier);
    }

    private final AllGather allGather;

    private class AllGather implements Runnable {

        final CyclicBarrier barrier = new CyclicBarrier(size, this);

        private final Object sendbuf[] = new Object[size],
                recvbuf[][] = new Object[size][size];

        private final int[] length = new int[size];

        public void setSendRecv(Object send, Object[] recv, int rank) {
            sendbuf[rank] = send;
            recvbuf[rank] = recv;
            length[rank] = Array.getLength(send);
        }

        public void run() {
            for (int i = 0; i < size; ++i)
                for (int j = 0; j < size; ++j)
                    System
                            .arraycopy(sendbuf[i], 0, recvbuf[j][i], 0,
                                    length[i]);
        }
    }

    void allToAll(Object[] sendbuf, Object[] recvbuf, int rank) {
        allToAll.setSendRecv(sendbuf, recvbuf, rank);

        await(allToAll.barrier);
    }

    private final AllToAll allToAll;

    private class AllToAll implements Runnable {

        final CyclicBarrier barrier = new CyclicBarrier(size, this);

        private final Object[][] sendbuf = new Object[size][size],
                recvbuf = new Object[size][size];

        private final int[][] length = new int[size][size];

        public void setSendRecv(Object[] send, Object[] recv, int rank) {
            sendbuf[rank] = send;
            recvbuf[rank] = recv;

            for (int i = 0; i < size; ++i)
                length[rank][i] = Array.getLength(send[i]);
        }

        public void run() {
            for (int i = 0; i < size; ++i)
                for (int j = 0; j < size; ++j)
                    System.arraycopy(sendbuf[i][j], 0, recvbuf[j][i], 0,
                            length[i][j]);
        }
    }

    void reduce(Object sendbuf, Object recvbuf, Reduction op, int root, int rank) {
        reduce.sendbuf[rank] = sendbuf;
        if (rank == root) {
            reduce.op = op;
            reduce.recvbuf = recvbuf;
        }

        await(reduce.barrier);
    }

    private final Reduce reduce;

    private class Reduce implements Runnable {

        CyclicBarrier barrier = new CyclicBarrier(size, this);

        Reduction op;

        Object sendbuf[] = new Object[size], recvbuf;

        public void run() {
            op.init(recvbuf);

            for (int i = 0; i < size; ++i)
                op.op(recvbuf, sendbuf[i]);
        }
    }

    void allReduce(Object sendbuf, Object recvbuf, Reduction op, int rank) {
        reduce(sendbuf, recvbuf, op, 0, rank);
        broadcast(recvbuf, 0, rank);
    }
}
