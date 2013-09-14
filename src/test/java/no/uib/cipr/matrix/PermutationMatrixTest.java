package no.uib.cipr.matrix;

import junit.framework.TestCase;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

/**
 * @author Sam Halliday
 */
public class PermutationMatrixTest extends TestCase {

  Matrix m = new DenseMatrix(new double[][]{
      {2, -1, -2},
      {-4, 6, 3},
      {-4, -2, 8}
  });
  DenseMatrix e = new DenseMatrix(new double[][]{
      {-4, 6, 3},
      {-4, -2, 8},
      {2, -1, -2}
  });
  DenseMatrix eI = new DenseMatrix(new double[][]{
      {-4, -2, 8},
      {2, -1, -2},
      {-4, 6, 3}
  });
  int[] piv = new int[]{2, 3, 3};

  public void testMultiply() {
    Matrix p = PermutationMatrix.fromPartialPivots(piv);
    Matrix c = p.mult(m, new CompRowMatrix(new DenseMatrix(m.numRows(), m.numColumns())));
    MatrixTestAbstract.assertEquals(e, c);
  }

  public void testMultiplyTrans() {
    Matrix p = PermutationMatrix.fromPartialPivots(piv);
    Matrix c = p.transAmult(m, new CompRowMatrix(new DenseMatrix(m.numRows(), m.numColumns())));
    MatrixTestAbstract.assertEquals(eI, c);
  }

  public void testMultiplyDense() {
    Matrix p = PermutationMatrix.fromPartialPivots(piv);
    Matrix c = p.mult(m, new DenseMatrix(m.numRows(), m.numColumns()));
    MatrixTestAbstract.assertEquals(e, c);
  }

  public void testMultiplyTransDense() {
    Matrix p = PermutationMatrix.fromPartialPivots(piv);
    Matrix c = p.transAmult(m, new DenseMatrix(m.numRows(), m.numColumns()));
    MatrixTestAbstract.assertEquals(eI, c);
  }

  public void testMultiplyCompRow() {
    Matrix p = PermutationMatrix.fromPartialPivots(piv);
    Matrix c = p.mult(new CompRowMatrix(m), new CompRowMatrix(new DenseMatrix(m.numRows(), m.numColumns())));
    MatrixTestAbstract.assertEquals(e, c);
  }

  public void testMultiplyTransCompRow() {
    Matrix p = PermutationMatrix.fromPartialPivots(piv);
    Matrix c = p.transAmult(new CompRowMatrix(m), new CompRowMatrix(new DenseMatrix(m.numRows(), m.numColumns())));
    MatrixTestAbstract.assertEquals(eI, c);
  }

}
