package no.uib.cipr.matrix.sparse;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;
import junit.framework.Assert;
import lombok.extern.java.Log;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
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

    for (MatrixEntry e : A) {
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

    LinkedSparseMatrix.Node node = ((LinkedSparseMatrix) A).links.headRow;
    while (node != null) {
//      log.info(node.toString());
      Assert.assertEquals(Ad[node.row][node.col], node.val);
      node = node.rowTail;
    }
    node = ((LinkedSparseMatrix) A).links.headCol;
    while (node != null) {
//      log.info(node.toString());
      Assert.assertEquals(Ad[node.row][node.col], node.val);
      node = node.colTail;
    }
  }

  public void ignoredTimedMult() {
    Stopwatch watch = Stopwatch.createUnstarted();
    DenseMatrix dense = new DenseMatrix(1000, 1000);
    int[][] nz = Utilities.getRowPattern(dense.numRows(), dense.numColumns(), 100);
    Utilities.rowPopulate(dense, nz);
    log.info("created matrices");
    Matrix sparse = new LinkedSparseMatrix(dense.numRows(), dense.numColumns());
    sparse.set(dense);

    for (Matrix m : Lists.newArrayList(dense, sparse)) {
      log.info("starting " + m.getClass());
      Matrix t = new DenseMatrix(m);
      t.transpose();
      Matrix o = new DenseMatrix(dense.numRows(), dense.numColumns());
      log.info("warming up " + m.getClass() + " " + o.getClass());
      for (int i = 0; i < 10; i++)
        m.mult(t, o);
      log.info("starting " + m.getClass() + " " + o.getClass());
      watch.start();
      for (int i = 0; i < 100; i++)
        m.mult(t, o);
      watch.stop();
      log.info(m.getClass() + " " + o.getClass() + " " + watch);
    }
  }

  public void ignoredTimedTransMult() {
    Stopwatch watch = Stopwatch.createUnstarted();
    DenseMatrix dense = new DenseMatrix(1000, 1000);
    int[][] nz = Utilities.getRowPattern(dense.numRows(), dense.numColumns(), 100);
    Utilities.rowPopulate(dense, nz);
    log.info("created matrices");
    Matrix sparse = new LinkedSparseMatrix(dense.numRows(), dense.numColumns());
    sparse.set(dense);

    for (Matrix m : Lists.newArrayList(dense, sparse)) {
      log.info("starting " + m.getClass());
      Matrix t = new DenseMatrix(m);
      Matrix o = new DenseMatrix(dense.numRows(), dense.numColumns());
      log.info("warming up " + m.getClass() + " " + o.getClass());
      for (int i = 0; i < 10; i++)
        m.transAmult(t, o);
      log.info("starting " + m.getClass() + " " + o.getClass());
      watch.start();
      for (int i = 0; i < 100; i++)
        m.transAmult(t, o);
      watch.stop();
      log.info(m.getClass() + " " + o.getClass() + " " + watch);
    }
  }

  @Override
  public void testIteratorSet() {
  }

  @Override
  public void testIteratorSetGet() {
  }
}
