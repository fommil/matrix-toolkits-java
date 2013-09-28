package no.uib.cipr.matrix.sparse;

import au.com.bytecode.opencsv.CSVWriter;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;
import junit.framework.Assert;
import lombok.Cleanup;
import lombok.extern.java.Log;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Utilities;

import java.io.File;
import java.io.FileWriter;
import java.util.concurrent.TimeUnit;

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

    LinkedSparseMatrix.Node head = ((LinkedSparseMatrix) A).links.head;
    LinkedSparseMatrix.Node node = head;
    while (node != null) {
//      log.info(node.toString());
      Assert.assertEquals(Ad[node.row][node.col], node.val);
      node = node.rowTail;
    }
    node = head;
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


  /**
   * Does a naive perf test against DenseMatrix, outputting CSV
   * that we plot in R.
   * <p/>
   * Generate n x n matrix with m entries, on left, and
   * an n x n dense matrix with m entries on the right,
   * using the same population algo. Then we multiply
   * them and output into a dense matrix. We compare
   * dense vs linked sparse from the left, and also
   * look at memory usage. We repeat for 10 different
   * values of m (10,000 to 100,000), and in both cases
   * vary m from 1,000 to 10,000. This requires about
   * 8GB heap to be on the safe side.
   * <p/>
   * -Xms8g -Xmx8g -Djava.util.logging.config.file=logging.properties
   *
   * @param args
   */
  public static void main(String[] args) throws Exception {
    File file = new File("LinkedSparseMatrixPerf.csv");
    log.info("writing to " + file);
    @Cleanup CSVWriter csv = new CSVWriter(new FileWriter(file));

    for (int r = 0; r < 10; r++) {
      for (int m = 10000; m <= 100000; m = m + 10000) {
        for (int n = 1000; n <= 10000; n = n + 1000) {
          int[][] patternA = Utilities.getRowPattern(n, n, m / n);
          DenseMatrix origA = new DenseMatrix(n, n);
          Utilities.rowPopulate(origA, patternA);
          int[][] patternB = Utilities.getRowPattern(n, n, m / n);
          DenseMatrix origB = new DenseMatrix(n, n);
          Utilities.rowPopulate(origB, patternB);
          // to be fair, we reuse the same matrix values

          long denseMem, denseInitTime, denseMultTime, sparseMem, sparseInitTime, sparseMultTime;

          Stopwatch timer = Stopwatch.createUnstarted();
          {
            timer.reset();
            timer.start();
            DenseMatrix A = new DenseMatrix(origA);
            timer.stop();
            // all attempts to measure memory usage failed
            denseMem = n * n * 8;
            denseInitTime = timer.elapsed(TimeUnit.NANOSECONDS);
            timer.reset();

            DenseMatrix B = origB.copy();
            DenseMatrix C = new DenseMatrix(n, n);
            timer.start();
            A.mult(B, C);
            timer.stop();
            denseMultTime = timer.elapsed(TimeUnit.NANOSECONDS);
          }
          {
            timer.reset();
            timer.start();
            LinkedSparseMatrix A = new LinkedSparseMatrix(origA);
            timer.stop();
            // using compressedooms
            sparseMem = m * 28 + 16 * n;
            sparseInitTime = timer.elapsed(TimeUnit.NANOSECONDS);
            timer.reset();

            DenseMatrix B = origB.copy();
            DenseMatrix C = new DenseMatrix(n, n);
            timer.start();
            A.mult(B, C);
            timer.stop();
            sparseMultTime = timer.elapsed(TimeUnit.NANOSECONDS);
          }

          String[] line = new String[]{
              Integer.toString(n), Integer.toString(m),
              Long.toString(denseMem), Long.toString(denseInitTime), Long.toString(denseMultTime),
              Long.toString(sparseMem), Long.toString(sparseInitTime), Long.toString(sparseMultTime)
          };
          log.info(java.util.Arrays.toString(line));
          csv.writeNext(line);

          // these are to keep lots of refs above alive from GC
          log.finest(origA.numRows() + " " + origB.numColumns() + " " + patternA.length + " " + patternB.length);
        }
      }
    }
  }


}
