package no.uib.cipr.matrix;

import java.util.Random;

/**
 * @author Sam Halliday
 */
public class DenseVectorSubTest extends VectorTestAbstract {

  public DenseVectorSubTest(String arg0) {
    super(arg0);
  }

  @Override
  protected void createPrimary() throws Exception {
    int n = Utilities.getInt(1, max);
    DenseVector wrapped = new DenseVector(n * 10);
    int offset = new Random().nextInt(n * 9);
    x = new DenseVectorSub(wrapped, offset, n);
    xd = Utilities.populate(x);
  }

}
