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
public class ArpackSymTest extends TestCase {

  public void testRandomEigensystem() throws NotConvergedException {
    for (int i = 100; i <= 500; i = i + 100) {
      UpperSymmDenseMatrix matrix = new UpperSymmDenseMatrix(i);
      Utilities.upperPopulateGauss(matrix);
      Map<Double, DenseVector> evd = evdSolve(matrix);

      ArpackSym solver = new ArpackSym(matrix);
      int todo = i / 10;

      Map<Double, DenseVectorSub> results = solver.solve(todo, ArpackSym.Ritz.LA);
      Assert.assertEquals(todo, results.size());
      for (Map.Entry<Double, DenseVectorSub> e : results.entrySet()) {
        // exact match of eigenvector / eigenvalue is not important for random matrices
        // as the eigenvectors should always be the Euclidean directions
        boolean value = false, vector = false;
        for (Map.Entry<Double, DenseVector> evdEntry : evd.entrySet()) {
          if (e.getKey() - evdEntry.getKey() < 0.0001)
            value = true;
          if (e.getValue().dot(evdEntry.getValue()) < 0.9999)
            vector = true;
        }
        if (!value)
          fail("no matching eigenvalue for " + e.getKey() + " : " + evd.keySet());
        if (!vector)
          fail("no matching eigenvector for " + e.getKey() + " : " + evd.keySet());
      }
    }
  }

  private Map<Double, DenseVector> evdSolve(UpperSymmDenseMatrix matrix) throws NotConvergedException {
    SymmDenseEVD evd = new SymmDenseEVD(matrix.numColumns(), true);
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

  public static void main(String[] args) throws Exception {
    File file = new File("A.txt");
    Matrix A = new LinkedSparseMatrix(new MatrixVectorReader(new FileReader(file)));
    Matrix AtA = A.transAmult(A, new LinkedSparseMatrix(A.numColumns(), A.numColumns()));

    ArpackSym solver = new ArpackSym(AtA);
    Map<Double, DenseVectorSub> results = solver.solve(A.numRows() / 10, ArpackSym.Ritz.LA);
    for (Map.Entry<Double, DenseVectorSub> result : results.entrySet()) {
      log.info(result.getKey().toString());
    }
  }

}
