package no.uib.cipr.matrix.sparse;

import com.github.fommil.netlib.ARPACK;
import lombok.RequiredArgsConstructor;
import lombok.extern.java.Log;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import org.netlib.util.doubleW;
import org.netlib.util.intW;

import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;

/**
 * Uses ARPACK to partially solve square symmetric eigensystems
 * (ARPACK is designed to compute a subset of eigenvalues/eigenvectors).
 *
 * @author Sam Halliday
 */
@RequiredArgsConstructor
@Log
public class ArpackSqSym {

  private final ARPACK arpack = ARPACK.getInstance();

  private static final double TOL = 0;

  private final Matrix matrix;

  /**
   * Solve the eigensystem for the number of eigenvalues requested.
   *
   * @param eigenvalues
   * @return a map from eigenvalues to corresponding eigenvectors.
   */
  public Map<Double, DenseVector> solve(int eigenvalues) {
    if (eigenvalues <= 0 || eigenvalues >= matrix.numColumns())
      throw new IllegalArgumentException(eigenvalues + " >= " + matrix.numColumns());

    int n = matrix.numRows();
    intW nev = new intW(eigenvalues);

    int ncv = Math.min(2 * eigenvalues, n);
    String bmat = "I";
    String which = "LA"; // prefer algebraically larger eigenvalues
    doubleW tol = new doubleW(TOL);
    intW info = new intW(0);
    int[] iparam = new int[11];
    iparam[0] = 1;
    iparam[2] = 300;
    iparam[6] = 1;
    intW ido = new intW(0);

    // used for initial residual (if info != 0)
    // and eventually the output residual
    double[] resid = new double[n];
    // Lanczos basis vectors
    double[] v = new double[n * ncv];
    // Arnoldi reverse communication
    double[] workd = new double[3 * n];
    // private work array
    double[] workl = new double[ncv * (ncv + 8)];
    int[] ipntr = new int[11];

    int i = 0;
    while (true) {
      i++;
      arpack.dsaupd(ido, bmat, n, which, nev.val, tol, resid, ncv, v, n, iparam, ipntr, workd, workl, workl.length, info);
      if (ido.val != -1 && ido.val != 1) break;
      // could be refactored to handle the other types of mode
      av(workd, ipntr[0] - 1, ipntr[1] - 1);
    }

    log.fine(i + " iterations for " + n);

    if (info.val < 0) throw new IllegalStateException("info = " + info.val);

    double[] d = new double[nev.val];
    boolean[] select = new boolean[ncv];
    double[] z = java.util.Arrays.copyOfRange(v, 0, nev.val * n);
    arpack.dseupd(true, "A", select, d, z, n, 0, bmat, n, which, nev, TOL, resid, ncv, v, n, iparam, ipntr, workd, workl, workl.length, info);
    if (info.val != 0) throw new IllegalStateException("info = " + info.val);

    int computed = iparam[4];
    log.fine("computed " + computed + " eigenvalues");

    Map<Double, DenseVector> solution = new TreeMap<Double, DenseVector>(new Comparator<Double>() {
      @Override
      public int compare(Double o1, Double o2) {
        // highest first
        return Double.compare(o2,o1);
      }
    });
    for (i = 0; i < computed ; i++) {
      double eigenvalue = d[i];
      double[] eigenvector = java.util.Arrays.copyOfRange(z, i * n, i * n + n);
      DenseVector vector = new DenseVector(eigenvector);
      solution.put(eigenvalue, vector);
    }

    return solution;
  }

  private void av(double[] work, int input_offset, int output_offset) {
    // RFE: https://github.com/fommil/matrix-toolkits-java/issues/35
    Vector x = new DenseVector(java.util.Arrays.copyOfRange(work, input_offset, input_offset + matrix.numColumns()));
    Vector y = new DenseVector(matrix.numColumns());
    matrix.mult(x, y);
    for (int i = 0; i < matrix.numColumns(); i++) {
      work[i + output_offset] = y.get(i);
    }
  }
}
