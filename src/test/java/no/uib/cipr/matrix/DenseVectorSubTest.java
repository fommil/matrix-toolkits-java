package no.uib.cipr.matrix;

import java.util.Random;

/**
 * @author Sam Halliday
 */
public class DenseVectorSubTest extends VectorTestAbstract {

    @Override
    protected void createPrimary() throws Exception {
        int n = Utilities.getInt(1, max);
        DenseVector wrapped = new DenseVector(n * 10);
        int offset = new Random().nextInt(n * 9);
        x = new DenseVectorSub(wrapped, offset, n);
        xd = Utilities.populate(x);
    }

}
