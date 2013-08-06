package no.uib.cipr.matrix.distributed;

import java.util.Random;

import junit.framework.TestCase;
import no.uib.cipr.matrix.distributed.CollectiveCommunications;
import no.uib.cipr.matrix.distributed.Communicator;

/**
 * Tests point-to-point communications of seven basic types
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
public class PointToPointTest extends TestCase {

    CollectiveCommunications coll;

    byte[] bsend, brecv;

    char[] csend, crecv;

    short[] ssend, srecv;

    int[] isend, irecv;

    long[] lsend, lrecv;

    float[] fsend, frecv;

    double[] dsend, drecv;

    int length;

    @Override
    protected void setUp() throws Exception {
        coll = new CollectiveCommunications(2);

        Random r = new Random();
        length = r.nextInt(100);

        bsend = new byte[length];
        brecv = new byte[length];

        csend = new char[length];
        crecv = new char[length];

        ssend = new short[length];
        srecv = new short[length];

        isend = new int[length];
        irecv = new int[length];

        lsend = new long[length];
        lrecv = new long[length];

        fsend = new float[length];
        frecv = new float[length];

        dsend = new double[length];
        drecv = new double[length];

        for (int i = 0; i < length; ++i) {
            dsend[i] = r.nextDouble();
            fsend[i] = r.nextFloat();
            lsend[i] = r.nextLong();
            isend[i] = r.nextInt();
            ssend[i] = (short) r.nextInt();
            csend[i] = (char) r.nextInt();
            bsend[i] = (byte) r.nextInt();
        }
    }

    public void testByteSendRecv() throws Exception {
        runSendRecv(bsend, brecv);

        for (int i = 0; i < length; ++i)
            assertEquals(bsend[i], brecv[i]);
    }

    public void testCharSendRecv() throws Exception {
        runSendRecv(csend, crecv);

        for (int i = 0; i < length; ++i)
            assertEquals(csend[i], crecv[i]);
    }

    public void testShortSendRecv() throws Exception {
        runSendRecv(ssend, srecv);

        for (int i = 0; i < length; ++i)
            assertEquals(ssend[i], srecv[i]);
    }

    public void testIntSendRecv() throws Exception {
        runSendRecv(isend, irecv);

        for (int i = 0; i < length; ++i)
            assertEquals(isend[i], irecv[i]);
    }

    public void testLongSendRecv() throws Exception {
        runSendRecv(lsend, lrecv);

        for (int i = 0; i < length; ++i)
            assertEquals(lsend[i], lrecv[i]);
    }

    public void testFloatSendRecv() throws Exception {
        runSendRecv(fsend, frecv);

        for (int i = 0; i < length; ++i)
            assertEquals(fsend[i], frecv[i], 1e-10);
    }

    public void testDoubleSendRecv() throws Exception {
        runSendRecv(dsend, drecv);

        for (int i = 0; i < length; ++i)
            assertEquals(dsend[i], drecv[i], 1e-10);
    }

    private void runSendRecv(Object send, Object recv)
            throws InterruptedException {
        Thread sender = createSender(send, 0, 1);
        Thread receiver = createReceiver(recv, 1, 0);

        sender.start();
        receiver.start();

        sender.join();
        receiver.join();
    }

    private Thread createSender(final Object send, final int rank,
            final int peer) {
        return new Thread(new Runnable() {
            public void run() {
                Communicator c = coll.createCommunicator(rank);
                c.send(send, peer);
            }
        });
    }

    private Thread createReceiver(final Object recv, final int rank,
            final int peer) {
        return new Thread(new Runnable() {
            public void run() {
                Communicator c = coll.createCommunicator(rank);
                c.recv(recv, peer);
            }
        });
    }
}
