package no.uib.cipr.matrix;

import junit.framework.Assert;
import junit.framework.TestCase;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

/**
 * @author Sam Halliday
 */
public class PermutationMatrixTest extends TestCase {

  public void testMultiply() {
    Matrix m = new CompRowMatrix(new DenseMatrix(new double[][]{
        {2, -1, -2},
        {-4, 6, 3},
        {-4, -2, 8}
    }));
    DenseMatrix e = new DenseMatrix(new double[][]{
        {-4, 6, 3},
        {-4, -2, 8},
        {2, -1, -2}
    });

    int[] piv = new int[]{2, 3, 3};

    Matrix p = PermutationMatrix.fromPartialPivots(piv);
    Matrix c = p.mult(m, new DenseMatrix(m.numRows(), m.numColumns()));

    for (int i = 0; i < m.numRows(); i++) {
      for (int j = 0; j < m.numColumns(); j++) {
        Assert.assertEquals(e.get(i, j), c.get(i, j));
      }
    }
  }

  public void testMultiplyDense() {
    DenseMatrix m = new DenseMatrix(new double[][]{
        {2, -1, -2},
        {-4, 6, 3},
        {-4, -2, 8}
    });
    DenseMatrix e = new DenseMatrix(new double[][]{
        {-4, 6, 3},
        {-4, -2, 8},
        {2, -1, -2}
    });

    int[] piv = new int[]{2, 3, 3};

    Matrix p = PermutationMatrix.fromPartialPivots(piv);
    Matrix c = p.mult(m, new DenseMatrix(m.numRows(), m.numColumns()));

    for (int i = 0; i < m.numRows; i++) {
      for (int j = 0; j < m.numColumns; j++) {
        Assert.assertEquals(e.get(i, j), c.get(i, j));
      }
    }
  }

}
