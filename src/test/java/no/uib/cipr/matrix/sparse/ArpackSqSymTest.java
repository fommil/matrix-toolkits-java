package no.uib.cipr.matrix.sparse;

import junit.framework.Assert;
import junit.framework.TestCase;
import lombok.extern.java.Log;
import no.uib.cipr.matrix.*;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Sam Halliday
 */
@Log
public class ArpackSqSymTest extends TestCase {


  public void testEigensystem() throws NotConvergedException {
    for (int i = 10; i < 1000; i = i + 10) {
      DenseMatrix matrix = new DenseMatrix(i, i);
      Utilities.populate(matrix);

      Map<Double, DenseVector> evd = evdSolve(matrix);

      ArpackSqSym solver = new ArpackSqSym(matrix);
      Map<Double, DenseVector> results = solver.solve(i - 1);

      Assert.assertEquals(i - 1, results.size());
      result: for (Map.Entry<Double, DenseVector> e : results.entrySet()) {
        for (Map.Entry<Double, DenseVector> evdEntry : evd.entrySet()) {
          log.info("eigenvalue = " + e.getKey());
          if (e.getKey() - evdEntry.getKey() < 0.0001) {
            MatrixTestAbstract.assertEquals(evdEntry.getValue(), e.getValue());
            continue result;
          }
          if (e.getKey() + evdEntry.getKey() < 0.0001) {
            DenseVector alt = new DenseVector(evdEntry.getValue());
            alt.scale(-1);
            MatrixTestAbstract.assertEquals(alt, e.getValue());
            continue result;
          }
          fail("no matching eigenvalues: " + evd.keySet());
        }
      }
    }
  }

  private Map<Double, DenseVector> evdSolve(DenseMatrix matrix) throws NotConvergedException {
    EVD evd = new EVD(matrix.numColumns());
    evd.factor(matrix);
    Map<Double, DenseVector> results = new HashMap<Double, DenseVector>();

    double[] eigenvalues = evd.getRealEigenvalues();
    Matrix eigenvectorMatrix = evd.getRightEigenvectors();
    for (int i = 0; i < eigenvalues.length; i++) {
      double eigenvalue = eigenvalues[i];
      DenseVector eigenvector = Matrices.getColumn(eigenvectorMatrix, i);
      double scale = eigenvector.norm(Vector.Norm.Two);
      log.info("scale = " + scale + ", would be " + (scale * eigenvalue));
      results.put(eigenvalue, eigenvector);
    }
    return results;
  }

}
