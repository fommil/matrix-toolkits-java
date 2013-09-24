package no.uib.cipr.matrix;

import java.util.Arrays;

/**
 * Wraps a DenseVector, allowing easy access to a sub array of
 * the original without taking copies.
 * <p/>
 * It should be possible to utilise BLAS / LAPACK in various
 * matrix classes. However, as it would be a mammoth task, it
 * will be done on an as-needed basis.
 *
 * @author Sam Halliday
 */
public class DenseVectorSub extends AbstractVector {

  private DenseVector wrapped;
  private int offset;

  public DenseVectorSub(DenseVector wrapped, int offset, int size) {
    super(size);
    if (offset + size > wrapped.size)
      throw new IllegalArgumentException(offset + "+" + size + ">" + wrapped.size);
    this.offset = offset;
    this.wrapped = wrapped;
  }

  @Override
  public double get(int index) {
    check(index);
    return wrapped.get(offset + index);
  }

  @Override
  public void set(int index, double value) {
    check(index);
    wrapped.set(offset + index, value);
  }

  @Override
  public DenseVector copy() {
    double[] data = Arrays.copyOfRange(wrapped.getData(), offset, offset + size);
    return new DenseVector(data, false);
  }
}
