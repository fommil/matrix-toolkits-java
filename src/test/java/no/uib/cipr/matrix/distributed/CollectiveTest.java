package no.uib.cipr.matrix.distributed;

import java.util.Random;

import no.uib.cipr.matrix.distributed.CollectiveCommunications;
import no.uib.cipr.matrix.distributed.Communicator;
import no.uib.cipr.matrix.distributed.Reductions;
import junit.framework.TestCase;

/**
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
public class CollectiveTest extends TestCase {

    CollectiveCommunications coll;

    Random r;

    @Override
    protected void setUp() throws Exception {
        r = new Random();
        int size = Math.max(1, r.nextInt(16));

        coll = new CollectiveCommunications(size);
    }

    public void testBarrier() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];

        for (int i = 0; i < t.length; ++i) {
            final int rank = i;
            t[i] = new Thread(new Runnable() {
                public void run() {
                    Communicator comm = coll.createCommunicator(rank);
                    comm.barrier();
                }
            });

            t[i].start();
        }

        for (int i = 0; i < t.length; ++i)
            t[i].join(10000);

        for (int i = 0; i < t.length; ++i)
            assertTrue(!t[i].isAlive());
    }

    public void testBroadcast() throws InterruptedException {
        final int[] send = new int[r.nextInt(100)];
        for (int i = 0; i < send.length; ++i)
            send[i] = r.nextInt();

        final int[][] recv = new int[coll.size()][send.length];
        recv[0] = send;

        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i) {
            final int rank = i;
            t[i] = new Thread(new Runnable() {
                public void run() {
                    Communicator comm = coll.createCommunicator(rank);
                    comm.broadcast(recv[rank], 0);
                }
            });

            t[i].start();
        }

        for (int i = 0; i < t.length; ++i)
            t[i].join();

        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < send.length; ++j)
                assertEquals(send[j], recv[i][j]);
    }

    public void testScatter() throws InterruptedException {
        int length = r.nextInt(100);
        final int[][] send = new int[coll.size()][length];
        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                send[i][j] = r.nextInt();

        final int[][] recv = new int[coll.size()][length];

        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i) {
            final int rank = i;
            t[i] = new Thread(new Runnable() {
                public void run() {
                    Communicator comm = coll.createCommunicator(rank);
                    comm.scatter(send, recv[rank], 0);
                }
            });

            t[i].start();
        }

        for (int i = 0; i < t.length; ++i)
            t[i].join();

        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                assertEquals(send[i][j], recv[i][j]);
    }

    public void testGather() throws InterruptedException {
        int length = r.nextInt(100);
        final int[][] send = new int[coll.size()][length];
        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                send[i][j] = r.nextInt();

        final int[][] recv = new int[coll.size()][length];

        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i) {
            final int rank = i;
            t[i] = new Thread(new Runnable() {
                public void run() {
                    Communicator comm = coll.createCommunicator(rank);
                    comm.gather(send[rank], recv, 0);
                }
            });

            t[i].start();
        }

        for (int i = 0; i < t.length; ++i)
            t[i].join();

        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                assertEquals(send[i][j], recv[i][j]);
    }

    public void testAllGather() throws InterruptedException {
        int length = r.nextInt(100);
        final int[][] send = new int[coll.size()][length];
        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                send[i][j] = r.nextInt();

        final int[][][] recv = new int[coll.size()][coll.size()][length];

        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i) {
            final int rank = i;
            t[i] = new Thread(new Runnable() {
                public void run() {
                    Communicator comm = coll.createCommunicator(rank);
                    comm.allGather(send[rank], recv[rank]);
                }
            });

            t[i].start();
        }

        for (int i = 0; i < t.length; ++i)
            t[i].join();

        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                assertEquals(send[i][j], recv[i][i][j]);
    }

    public void testAllToAll() throws InterruptedException {
        int length = r.nextInt(100);
        final int[][][] send = new int[coll.size()][coll.size()][length];
        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < coll.size(); ++j)
                for (int k = 0; k < length; ++k)
                    send[i][j][k] = r.nextInt();

        final int[][][] recv = new int[coll.size()][coll.size()][length];

        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i) {
            final int rank = i;
            t[i] = new Thread(new Runnable() {
                public void run() {
                    Communicator comm = coll.createCommunicator(rank);
                    comm.allToAll(send[rank], recv[rank]);
                }
            });

            t[i].start();
        }

        for (int i = 0; i < t.length; ++i)
            t[i].join();

        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < coll.size(); ++j)
                for (int k = 0; k < length; ++k)
                    assertEquals(send[i][j][k], recv[j][i][k]);
    }

    public void testReduce() throws InterruptedException {
        int length = r.nextInt(100);
        final int[][] send = new int[coll.size()][length];
        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                send[i][j] = r.nextInt(1000);

        final int[] recv = new int[length];

        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i) {
            final int rank = i;
            t[i] = new Thread(new Runnable() {
                public void run() {
                    Communicator comm = coll.createCommunicator(rank);
                    comm.reduce(send[rank], recv, Reductions.sum(), 0);
                }
            });

            t[i].start();
        }

        for (int i = 0; i < t.length; ++i)
            t[i].join();

        int[] mysum = new int[length];
        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                mysum[j] += send[i][j];

        for (int j = 0; j < length; ++j)
            assertEquals(mysum[j], recv[j]);
    }

    public void testAllReduce() throws InterruptedException {
        int length = r.nextInt(100);
        final int[][] send = new int[coll.size()][length];
        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                send[i][j] = r.nextInt(1000);

        final int[][] recv = new int[coll.size()][length];

        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i) {
            final int rank = i;
            t[i] = new Thread(new Runnable() {
                public void run() {
                    Communicator comm = coll.createCommunicator(rank);
                    comm.allReduce(send[rank], recv[rank], Reductions.sum());
                }
            });

            t[i].start();
        }

        for (int i = 0; i < t.length; ++i)
            t[i].join();

        int[] mysum = new int[length];
        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                mysum[j] += send[i][j];

        for (int i = 0; i < coll.size(); ++i)
            for (int j = 0; j < length; ++j)
                assertEquals(mysum[j], recv[i][j]);
    }
}
