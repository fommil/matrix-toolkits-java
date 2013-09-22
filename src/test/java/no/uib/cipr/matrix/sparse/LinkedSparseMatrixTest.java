package no.uib.cipr.matrix.sparse;

import junit.framework.Assert;
import lombok.extern.java.Log;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Utilities;

/**
 * @author Sam Halliday
 */
@Log
public class LinkedSparseMatrixTest extends SparseMatrixTestAbstract {
  public LinkedSparseMatrixTest(String arg0) {
    super(arg0);
  }

  @Override
  protected void createPrimary() throws Exception {
    int n = Utilities.getInt(1, max);
    int m = Utilities.getInt(1, max);
    int b = Utilities.getInt(Math.min(bmax, m));
    int[][] nz = Utilities.getRowPattern(n, m, b);
    A = new LinkedSparseMatrix(n, m);
    Ad = Utilities.rowPopulate(A, nz);

    for(MatrixEntry e : A) {
      int row = e.row();
      int col = e.column();
      double expect = Ad[row][col];
      Assert.assertEquals(expect, e.get(), 0);
      Assert.assertEquals(expect, A.get(row, col), 0);
    }

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        double expect = Ad[i][j];
        Assert.assertEquals(expect, A.get(i, j), 0);
      }
    }
  }

  @Override
  public void testIteratorSet() {
  }

  @Override
  public void testIteratorSetGet() {
  }

  @Override
  public void testCopy() {
  }
}
