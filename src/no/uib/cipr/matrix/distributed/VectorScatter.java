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

import java.util.concurrent.Future;

import no.uib.cipr.matrix.Vector;

/**
 * Vector scatter. Transfers values between distributed vectors and local
 * vectors, both ways (scatters and gathers). After starting a scatter or
 * gather, it must be ended using one of the two possible scatters/gathers
 * closing operations (set or add).
 * <p>
 * The operations in this class are not thread-safe.
 * </p>
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
class VectorScatter {

    /**
     * Communicator in use
     */
    private final Communicator comm;

    /**
     * For asynchronous communications
     */
    private final Future[] t;

    /**
     * Contains the entries recieved from other ranks. Must be used in
     * conjunction with the indices to assemble the local vector
     */
    private final double[][] recvD;

    /**
     * Contains the indices to put entries recieved from other ranks into the
     * local vector. Used during local vector assembly
     */
    private final int[][] recvI;

    /**
     * The indices of the local part of the global vector to send to the other
     * ranks
     */
    private final int[][] sendI;

    /**
     * Will contain the data of local part of the global vector which is to be
     * sent to the other ranks
     */
    private final double[][] sendD;

    /**
     * Who to send/receive to/from
     */
    private final boolean[] send, recv;

    /**
     * Constructor for VecScatter
     * 
     * @param comm
     *            Communicator to operate within
     * @param sendI
     *            The indices to send data from to each rank
     * @param recvI
     *            The indices to recieve data from each rank
     */
    public VectorScatter(Communicator comm, int[][] sendI, int[][] recvI) {
        this.comm = comm;
        this.sendI = sendI;
        this.recvI = recvI;

        if (sendI[comm.rank()].length != 0)
            throw new IllegalArgumentException("Illegal local send detected");
        if (recvI[comm.rank()].length != 0)
            throw new IllegalArgumentException("Illegal local recv detected");

        // Some preallocations
        t = new Future[2 * comm.size()];
        sendD = new double[comm.size()][];
        recvD = new double[comm.size()][];
        for (int i = 0; i < comm.size(); ++i) {
            sendD[i] = new double[sendI[i].length];
            recvD[i] = new double[recvI[i].length];
        }

        // For optimized communications
        send = new boolean[comm.size()];
        recv = new boolean[comm.size()];
        for (int i = 0; i < comm.size(); ++i) {
            send[i] = sendI[i].length > 0;
            recv[i] = recvI[i].length > 0;
        }
    }

    /**
     * Starts a scatter from x into y. It must be ended using an appropriate
     * operation
     */
    public void startScatter(DistVector x, Vector y) {
        if (x.size() != y.size())
            throw new IllegalArgumentException(
                    "Vectors must be of the same global size");

        // Assemble the arrays to send
        for (int i = 0; i < comm.size(); ++i)
            if (i != comm.rank())
                for (int j = 0; j < sendI[i].length; ++j)
                    sendD[i][j] = x.get(sendI[i][j]);

        // Start the sends and recieves
        for (int i = 0; i < comm.size(); ++i)
            if (i != comm.rank()) {
                if (send[i])
                    t[i] = comm.isend(sendD[i], i);
                if (recv[i])
                    t[i + comm.size()] = comm.irecv(recvD[i], i);
            }
    }

    /**
     * Finishes the scatter by inserting the values over those already in y
     */
    public void endSetScatter(DistVector x, Vector y) {
        if (x.size() != y.size())
            throw new IllegalArgumentException(
                    "Vectors must be of the same global size");

        // Finish pending communications
        comm.await(t);

        // Assemble the local vector by insertion
        for (int i = 0; i < comm.size(); ++i)
            if (i != comm.rank())
                for (int j = 0; j < recvI[i].length; ++j)
                    y.set(recvI[i][j], recvD[i][j]);
    }

    /**
     * Finished the scatter by adding the values to those already in y
     */
    public void endAddScatter(DistVector x, Vector y) {
        if (x.size() != y.size())
            throw new IllegalArgumentException(
                    "Vectors must be of the same global size");

        // Finish pending communications
        comm.await(t);

        // Assemble the local vector by addition
        for (int i = 0; i < comm.size(); ++i)
            if (i != comm.rank())
                for (int j = 0; j < recvI[i].length; ++j)
                    y.add(recvI[i][j], recvD[i][j]);
    }

    /**
     * Starts a gather from x into y. It must be ended using an apropriate
     * operation
     */
    public void startGather(Vector x, DistVector y) {
        if (x.size() != y.size())
            throw new IllegalArgumentException(
                    "Vectors must be of the same global size");

        // Assemble the arrays to send
        for (int i = 0; i < comm.size(); ++i)
            if (i != comm.rank())
                for (int j = 0; j < recvI[i].length; ++j)
                    recvD[i][j] = x.get(recvI[i][j]);

        // Start the sends and recieves
        for (int i = 0; i < comm.size(); ++i)
            if (i != comm.rank()) {
                if (send[i])
                    t[i] = comm.isend(recvD[i], i);
                if (recv[i])
                    t[i + comm.size()] = comm.irecv(sendD[i], i);
            }
    }

    /**
     * Finishes the gather by inserting the values over those already in y
     */
    public void endSetGather(Vector x, DistVector y) {
        if (x.size() != y.size())
            throw new IllegalArgumentException(
                    "Vectors must be of the same global size");

        // Finish pending communications
        comm.await(t);

        // Assemble the global vector by insertion
        for (int i = 0; i < comm.size(); ++i)
            if (i != comm.rank())
                for (int j = 0; j < sendI[i].length; ++j)
                    y.set(sendI[i][j], sendD[i][j]);
    }

    /**
     * Finishes the gather by adding the values to those already in y
     */
    public void endAddGather(Vector x, DistVector y) {
        if (x.size() != y.size())
            throw new IllegalArgumentException(
                    "Vectors must be of the same global size");

        // Finish pending communications
        comm.await(t);

        // Assemble the global vector by insertion
        for (int i = 0; i < comm.size(); ++i)
            if (i != comm.rank())
                for (int j = 0; j < sendI[i].length; ++j)
                    y.add(sendI[i][j], sendD[i][j]);
    }

}
