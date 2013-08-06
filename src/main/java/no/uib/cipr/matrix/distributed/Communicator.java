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
import java.util.List;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.Exchanger;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;

/**
 * Inter-thread communications. Supports point-to-point communications using
 * barriers between the threads. Construct it using the
 * <code>CollectiveCommunications.createCommunicator</code> method.
 * <p>
 * All <code>Object</code>s which are sent and recieved are arrays (for
 * instance, <code>double[]</code> or <code>int[]</code>), and the types
 * must be compatible. It follows that <code>Object[]</code> is an array of
 * native arrays, such as <code>int[][]</code>.
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
public class Communicator {

    /**
     * My rank
     */
    private final int rank;

    /**
     * Executes asynchronous operations
     */
    private final ExecutorService executor;

    /**
     * Collective communications
     */
    private final CollectiveCommunications coll;

    /**
     * Input and output locks for irecv/isend
     */
    final Object[] in, out;

    /**
     * Barrier for communication with a peer
     */
    final SendRecv[] send, recv;

    static class SendRecv implements Runnable {

        CyclicBarrier barrier = new CyclicBarrier(2, this);

        Object send, recv;

        int sendOffset, recvOffset, length;

        public void run() {
            System.arraycopy(send, sendOffset, recv, recvOffset, length);
        }
    }

    /**
     * Sets up a communicator between the given number of threads
     */
    Communicator(int rank, final List<Exchanger<SendRecv>> ex,
            CollectiveCommunications coll) {

        this.rank = rank;
        this.coll = coll;

        if (rank < 0)
            throw new IllegalArgumentException("rank < 0");
        if (rank >= size())
            throw new IllegalArgumentException("rank >= size");

        // Create daemon threads for running async communications
        executor = Executors.newCachedThreadPool(new ThreadFactory() {
            public Thread newThread(Runnable r) {
                Thread t = new Thread(r);
                t.setDaemon(true);
                return t;
            }
        });

        in = new Object[size()];
        out = new Object[size()];
        for (int i = 0; i < size(); ++i) {
            in[i] = new Object();
            out[i] = new Object();
        }

        send = new SendRecv[size()];
        recv = new SendRecv[size()];

        // Create local sends, and get recvs from the peers
        for (int i = 0; i < size(); ++i)
            if (i != rank)
                try {
                    send[i] = new SendRecv();
                    recv[i] = ex.get(i).exchange(send[i]);
                } catch (InterruptedException e) {
                    throw new RuntimeException(e);
                }
    }

    /**
     * Rank of this thread in the collective
     */
    public int rank() {
        return rank;
    }

    /**
     * Size of the collective
     */
    public int size() {
        return coll.size();
    }

    /**
     * Gathers data from all tasks and distribute it to all.
     * <p>
     * Row corresponds to ranks, and columns corresponds to data.
     * <p>
     * Input: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>B1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>C1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>D1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * </table>
     * <p>
     * Output: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>B1</td>
     * <td>C1</td>
     * <td>D1</td>
     * </tr>
     * <tr>
     * <td>A1</td>
     * <td>B1</td>
     * <td>C1</td>
     * <td>D1</td>
     * </tr>
     * <tr>
     * <td>A1</td>
     * <td>B1</td>
     * <td>C1</td>
     * <td>D1</td>
     * </tr>
     * <tr>
     * <td>A1</td>
     * <td>B1</td>
     * <td>C1</td>
     * <td>D1</td>
     * </tr>
     * </table>
     */
    public void allGather(Object sendbuf, Object[] recvbuf) {
        coll.allGather(sendbuf, recvbuf, rank);
    }

    /**
     * Combines values from all processes and distribute the result back to all
     * processes.
     */
    public void allReduce(Object sendbuf, Object recvbuf, Reduction op) {
        coll.allReduce(sendbuf, recvbuf, op, rank);
    }

    /**
     * Sends data from all to all processes.
     * <p>
     * Row corresponds to ranks, and columns corresponds to data.
     * <p>
     * Input: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>A2</td>
     * <td>A3</td>
     * <td>A4</td>
     * </tr>
     * <tr>
     * <td>B1</td>
     * <td>B2</td>
     * <td>B3</td>
     * <td>B4</td>
     * </tr>
     * <tr>
     * <td>C1</td>
     * <td>C2</td>
     * <td>C3</td>
     * <td>C4</td>
     * </tr>
     * <tr>
     * <td>D1</td>
     * <td>D2</td>
     * <td>D3</td>
     * <td>D4</td>
     * </tr>
     * </table>
     * <p>
     * Output: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>B1</td>
     * <td>C1</td>
     * <td>D1</td>
     * </tr>
     * <tr>
     * <td>A2</td>
     * <td>B2</td>
     * <td>C2</td>
     * <td>D2</td>
     * </tr>
     * <tr>
     * <td>A3</td>
     * <td>B3</td>
     * <td>C3</td>
     * <td>D3</td>
     * </tr>
     * <tr>
     * <td>A4</td>
     * <td>B4</td>
     * <td>C4</td>
     * <td>D4</td>
     * </tr>
     * </table>
     */
    public void allToAll(Object[] sendbuf, Object[] recvbuf) {
        coll.allToAll(sendbuf, recvbuf, rank);
    }

    /**
     * Blocks until all process have reached this routine.
     */
    public void barrier() {
        coll.barrier();
    }

    /**
     * Broadcasts a message from the process with rank "root" to all other
     * processes of the group.
     * <p>
     * Row corresponds to ranks, and columns corresponds to data.
     * <p>
     * Input: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * </table>
     * <p>
     * Output: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>A1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>A1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>A1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * </table>
     */
    public void broadcast(Object buffer, int root) {
        coll.broadcast(buffer, root, rank);
    }

    /**
     * Gathers together values from a group of processes.
     * <p>
     * Row corresponds to ranks, and columns corresponds to data.
     * <p>
     * Input: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>A2</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>A3</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>A4</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * </table>
     * <p>
     * Output: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>A2</td>
     * <td>A3</td>
     * <td>A4</td>
     * </tr>
     * <tr>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * </table>
     */
    public void gather(Object sendbuf, Object[] recvbuf, int root) {
        coll.gather(sendbuf, recvbuf, root, rank);
    }

    /**
     * Reduces values on all processes to a single value
     */
    public void reduce(Object sendbuf, Object recvbuf, Reduction op, int root) {
        coll.reduce(sendbuf, recvbuf, op, root, rank);
    }

    /**
     * Sends data from one task to all other tasks in a group.
     * <p>
     * Row corresponds to ranks, and columns corresponds to data.
     * <p>
     * Input: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>A2</td>
     * <td>A3</td>
     * <td>A4</td>
     * </tr>
     * <tr>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * </table>
     * <p>
     * Output: <table border="1">
     * <tr>
     * <td>A1</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>A2</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>A3</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * <tr>
     * <td>A4</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * <td>&nbsp;</td>
     * </tr>
     * </table>
     */
    public void scatter(Object[] sendbuf, Object recvbuf, int root) {
        coll.scatter(sendbuf, recvbuf, root, rank);
    }

    /**
     * Sends data[offset:offset+length] to peer
     */
    public void send(Object data, int offset, int length, int peer) {
        checkArgs(data, offset, length, peer);

        send[peer].length = length;
        send[peer].sendOffset = offset;
        send[peer].send = data;

        CollectiveCommunications.await(send[peer].barrier);
    }

    /**
     * Receives data[offset:offset+length] from peer
     */
    public void recv(Object data, int offset, int length, int peer) {
        checkArgs(data, offset, length, peer);

        recv[peer].recvOffset = offset;
        recv[peer].recv = data;

        CollectiveCommunications.await(recv[peer].barrier);
    }

    public Future isend(final Object data, final int offset, final int length,
            final int peer) {
        return executor.submit(new Runnable() {
            public void run() {
                synchronized (out[peer]) {
                    send(data, offset, length, peer);
                }
            }
        });
    }

    public Future irecv(final Object data, final int offset, final int length,
            final int peer) {
        return executor.submit(new Runnable() {
            public void run() {
                synchronized (in[peer]) {
                    recv(data, offset, length, peer);
                }
            }
        });
    }

    public void send(Object data, int peer) {
        send(data, 0, Array.getLength(data), peer);
    }

    public void recv(Object data, int peer) {
        recv(data, 0, Array.getLength(data), peer);
    }

    public Future isend(Object data, int peer) {
        return isend(data, 0, Array.getLength(data), peer);
    }

    public Future irecv(Object data, int peer) {
        return irecv(data, 0, Array.getLength(data), peer);
    }

    /**
     * Waits for the given asynchronous operations to finish
     */
    public void await(Future[] future) {
        for (Future f : future)
            await(f);
    }

    /**
     * Waits for the given asynchronous operation to finish
     */
    public void await(Future f) {
        if (f == null)
            return;
        try {
            f.get();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private void checkArgs(Object data, int offset, int length, int peer) {
        if (peer == rank)
            throw new IllegalArgumentException("peer == rank");
        if (length + offset > Array.getLength(data))
            throw new IllegalArgumentException("Buffer underflow");
        if (peer < 0 || peer >= coll.size)
            throw new IllegalArgumentException("Invalid peer");
    }
}
