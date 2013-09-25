package no.uib.cipr.matrix.sparse;

import com.github.fommil.netlib.ARPACK;
import lombok.extern.java.Log;
import no.uib.cipr.matrix.*;
import org.netlib.util.doubleW;
import org.netlib.util.intW;

import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;

/**
 * Uses ARPACK to partially solve symmetric eigensystems
 * (ARPACK is designed to compute a subset of eigenvalues/eigenvectors).
 *
 * @author Sam Halliday
 */
@Log
public class ArpackSym {

  public enum Ritz {
    /**
     * compute the NEV largest (algebraic) eigenvalues.
     */
    LA,
    /**
     * compute the NEV smallest (algebraic) eigenvalues.
     */
    SA,
    /**
     * compute the NEV largest (in magnitude) eigenvalues.
     */
    LM,
    /**
     * compute the NEV smallest (in magnitude) eigenvalues.
     */
    SM,
    /**
     * compute NEV eigenvalues, half from each end of the spectrum
     */
    BE
  }

  private final ARPACK arpack = ARPACK.getInstance();

  private static final double TOL = 0.0001;

  private static final boolean EXPENSIVE_CHECKS = true;

  private final Matrix matrix;


  public ArpackSym(Matrix matrix) {
    if (!matrix.isSquare())
      throw new IllegalArgumentException("matrix must be square");
    if (EXPENSIVE_CHECKS)
      for (MatrixEntry entry : matrix) {
        if (entry.get() != matrix.get(entry.column(), entry.row()))
          throw new IllegalArgumentException("matrix must be symmetric");
      }
    this.matrix = matrix;
  }


  /**
   * Solve the eigensystem for the number of eigenvalues requested.
   * <p>
   * NOTE: The references to the eigenvectors will keep alive a reference to
   * a {@code nev * n} double array, so use the {@code copy()} method to free
   * it up if only a subset is required.
   *
   * @param eigenvalues
   * @param ritz        preference for solutions
   * @return a map from eigenvalues to corresponding eigenvectors.
   */
  public Map<Double, DenseVectorSub> solve(int eigenvalues, Ritz ritz) {
    if (eigenvalues <= 0)
      throw new IllegalArgumentException(eigenvalues + " <= 0");
    if (eigenvalues >= matrix.numColumns())
      throw new IllegalArgumentException(eigenvalues + " >= " + (matrix.numColumns()));

    int n = matrix.numRows();
    intW nev = new intW(eigenvalues);

    int ncv = Math.min(2 * eigenvalues, n);

    String bmat = "I";
    String which = ritz.name();
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
      if (ido.val == 99) break;
      if (ido.val != -1 && ido.val != 1) throw new IllegalStateException("ido = " + ido.val);
      // could be refactored to handle the other types of mode
      av(workd, ipntr[0] - 1, ipntr[1] - 1);
    }

    ArpackSym.log.fine(i + " iterations for " + n);

    if (info.val != 0) throw new IllegalStateException("info = " + info.val);

    double[] d = new double[nev.val];
    boolean[] select = new boolean[ncv];
    double[] z = java.util.Arrays.copyOfRange(v, 0, nev.val * n);

    arpack.dseupd(true, "A", select, d, z, n, 0, bmat, n, which, nev, TOL, resid, ncv, v, n, iparam, ipntr, workd, workl, workl.length, info);
    if (info.val != 0) throw new IllegalStateException("info = " + info.val);

    int computed = iparam[4];
    ArpackSym.log.fine("computed " + computed + " eigenvalues");

    Map<Double, DenseVectorSub> solution = new TreeMap<Double, DenseVectorSub>(new Comparator<Double>() {
      @Override
      public int compare(Double o1, Double o2) {
        // highest first
        return Double.compare(o2, o1);
      }
    });
    DenseVector eigenvectors = new DenseVector(z, false);
    for (i = 0; i < computed; i++) {
      double eigenvalue = d[i];
      DenseVectorSub eigenvector = new DenseVectorSub(eigenvectors, i * n, n);
      solution.put(eigenvalue, eigenvector);
    }

    return solution;
  }

  private void av(double[] work, int input_offset, int output_offset) {
    DenseVector w = new DenseVector(work, false);
    Vector x = new DenseVectorSub(w, input_offset, matrix.numColumns());
    Vector y = new DenseVectorSub(w, output_offset, matrix.numColumns());
    matrix.mult(x, y);
  }
}
