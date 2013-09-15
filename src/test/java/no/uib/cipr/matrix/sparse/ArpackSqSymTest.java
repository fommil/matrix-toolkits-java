package no.uib.cipr.matrix.sparse;

import junit.framework.Assert;
import junit.framework.TestCase;
import lombok.extern.java.Log;
import no.uib.cipr.matrix.*;
import no.uib.cipr.matrix.io.MatrixVectorReader;

import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Sam Halliday
 */
@Log
public class ArpackSqSymTest extends TestCase {

  public void testRandomEigensystem() throws NotConvergedException {
    for (int i = 100; i <= 1000; i = i + 100) {
      UpperSymmDenseMatrix matrix = new UpperSymmDenseMatrix(i);
      Utilities.upperPopulateGauss(matrix);

      Map<Double, DenseVector> evd = evdSolve(matrix);

      ArpackSqSym solver = new ArpackSqSym(matrix);
      int todo = i / 10;

      Map<Double, DenseVector> results = solver.solve(todo);
      Assert.assertEquals(todo, results.size());
      for (Map.Entry<Double, DenseVector> e : results.entrySet()) {
        // exact match of eigenvector / eigenvalue is not important for random matrices
        // as the eigenvectors should always be the Euclidean directions
        boolean value = false, vector = false;
        for (Map.Entry<Double, DenseVector> evdEntry : evd.entrySet()) {
          if (e.getKey() - Math.abs(evdEntry.getKey()) < 0.0001)
            value = true;
          if (e.getValue().dot(evdEntry.getValue()) < 0.9999)
            vector = true;
        }
        if (!value)
          fail("no matching eigenvalue for " + e.getKey() +" : " + evd.keySet());
        if (!vector)
          fail("no matching eigenvector for " + e.getKey() +" : " + evd.keySet());
      }
    }
  }

  private Map<Double, DenseVector> evdSolve(UpperSymmDenseMatrix matrix) throws NotConvergedException {
    SymmDenseEVD evd= new SymmDenseEVD(matrix.numColumns(), true);
    evd.factor(matrix);
    Map<Double, DenseVector> results = new HashMap<Double, DenseVector>();

    double[] eigenvalues = evd.getEigenvalues();
    Matrix eigenvectorMatrix = evd.getEigenvectors();
    for (int i = 0; i < eigenvalues.length; i++) {
      double eigenvalue = eigenvalues[i];
      DenseVector eigenvector = Matrices.getColumn(eigenvectorMatrix, i);
      results.put(eigenvalue, eigenvector);
    }
    return results;
  }

  public static void main(String [] args) throws Exception {
    File file = new File("A.txt");
    Matrix m = new DenseMatrix(new MatrixVectorReader(new FileReader(file)));

    ArpackSqSym solver = new ArpackSqSym(m);
    Map<Double, DenseVector> results = solver.solve(m.numColumns() / 10);
    for (Map.Entry<Double, DenseVector> result: results.entrySet()) {
      log.info(result.getKey().toString());
    }
  }

}
