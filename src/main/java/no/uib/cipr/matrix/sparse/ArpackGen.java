package no.uib.cipr.matrix.sparse;

import com.github.fommil.netlib.ARPACK;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import org.netlib.util.doubleW;
import org.netlib.util.intW;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;

/**
 * Uses ARPACK to partially solve general eigensystems (ARPACK is designed to
 * compute a subset of eigenvalues/eigenvectors).
 *
 * @author Andreas Solti adopted from the class {@link no.uib.cipr.matrix.sparse.ArpackSym} by Sam Halliday
 */
public class ArpackGen {

    private static final java.util.logging.Logger log = java.util.logging.Logger.getLogger(ArpackGen.class.getName());
    private boolean computeOnlyEigenvalues;

    public void setComputeOnlyEigenvalues(boolean computeOnlyEigenvalues) {
        this.computeOnlyEigenvalues = computeOnlyEigenvalues;
    }


    public enum Ritz {
        /**
         * compute the NEV largest (in magnitude) eigenvalues.
         */
        LM,
        /**
         * compute the NEV smallest (in magnitude) eigenvalues.
         */
        SM,
        /**
         * the NEV eigenvalues of largest real part.
         */
        LR,
        /**
         * the NEV eigenvalues of smallest real part.
         */
        SR,
        /**
         * want the NEV eigenvalues of largest imaginary part.
         */
        LI,
        /**
         *  want the NEV eigenvalues of smallest imaginary part.
         */
        SI
    }

    private final ARPACK arpack = ARPACK.getInstance();

    private static final double TOLERANCE = 0.00001;
    private static final int MAX_ITERATIONS = 10000;

    private final Matrix matrix;

    /**
     * If the matrix is symmetric, {@link no.uib.cipr.matrix.sparse.ArpackSym} is the better choice, as this property
     * can be exploited in computation and
     * @param matrix
     */
    public ArpackGen(Matrix matrix) {
        if (!matrix.isSquare())
            throw new IllegalArgumentException("matrix must be square");
        this.matrix = matrix;
        this.computeOnlyEigenvalues = true;
    }

    /**
     * Solve the eigensystem for the number of eigenvalues requested.
     * <p>
     * NOTE: The references to the eigenvectors will keep alive a reference to a
     * {@code nev * n} double array, so use the {@code copy()} method to free it
     * up if only a subset is required.
     *
     * @param eigenvalues
     * @param ritz
     *            preference for solutions
     * @return a map from eigenvalues to corresponding eigenvectors.
     */
    public Map<Double, DenseVectorSub> solve(int eigenvalues, Ritz ritz) {
        return solve(eigenvalues, ritz, TOLERANCE);
    }

    /**
     * Solve the eigensystem for the number of eigenvalues requested.
     * <p>
     * NOTE: The references to the eigenvectors will keep alive a reference to a
     * {@code nev * n} double array, so use the {@code copy()} method to free it
     * up if only a subset is required.
     *
     * @param eigenvalues
     * @param ritz
     *            preference for solutions
     * @return a map from eigenvalues to corresponding eigenvectors.
     */
    public Map<Double, DenseVectorSub> solve(int eigenvalues, Ritz ritz, double tolerance) {
        return solve(eigenvalues, ritz, tolerance, 0);
    }

    /**
     *
     * @param eigenvalues  NEV     Integer.  (INPUT)
    c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1. (Where N is the matrix dimension)
     * @param ritz
     * @param tolerance
     * @param ncvModification
     *    At present there is no a-priori analysis to guide the selection
     *    of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
     *    However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
     *    the same type are to be solved, one should experiment with increasing
     *    NCV while keeping NEV fixed for a given test problem.  This will
     *    usually decrease the required number of OP*x operations but it
     *    also increases the work and storage required to maintain the orthogonal
     *    basis vectors.  The optimal "cross-over" with respect to CPU time
     *    is problem dependent and must be determined empirically.
     *    See Chapter 8 of Reference 2 for further information.
     *
     * @return
     */
    public Map<Double, DenseVectorSub> solve(int eigenvalues, Ritz ritz, double tolerance, int ncvModification) {
        if (eigenvalues <= 0)
            throw new IllegalArgumentException(eigenvalues + " <= 0");
        if (eigenvalues > matrix.numColumns() - 2 )
            throw new IllegalArgumentException("NEV " + eigenvalues + " >= "
                    + (matrix.numColumns())+" - 2");

        int n = matrix.numRows();
        intW nev = new intW(eigenvalues);

        int ncv = Math.min(2 + eigenvalues + ncvModification, n); // Number of columns of the matrix V. NCV must satisfy the two inequalities 2 <= NCV-NEV and NCV <= N.

        String bmat = "I";
        String which = ritz.name();
        doubleW tol = new doubleW(tolerance);
        intW info = new intW(0);
        int[] iparam = new int[11];
        iparam[0] = 1;
        iparam[2] = Math.max(MAX_ITERATIONS, n); // maximum number of Arnoldi update operations allowed (default = 300?)
        iparam[4] = 1; // IPARAM(4) = NB: blocksize to be used in the recurrence. The code currently works only for NB = 1.
        iparam[6] = 1; // No longer referenced. Implicit restarting is ALWAYS used.
        iparam[7] = 1; // MODE On INPUT determines what type of eigenproblem is being solved. Must be 1,2,3,4; See under \Description of dnaupd  for the four modes available. (Mode 1:  A*x = lambda*x)
        intW ido = new intW(0);

        // used for initial residual (if info != 0)
        // and eventually the output residual
        double[] resid = new double[n];
        // Lanczos basis vectors
        double[] v = new double[n * ncv];
        // Arnoldi reverse communication
        double[] workd = new double[3 * n];
        // private work array
        double[] workl = new double[3*ncv * ncv + 6*ncv];
        int[] ipntr = new int[11];

        int i = 0;
        while (true) {
            i++;
            arpack.dnaupd(ido, bmat, n, which, nev.val, tol, resid, ncv, v, n,
                    iparam, ipntr, workd, workl, workl.length, info);
            if (ido.val == 99){
                if (info.val != 1) {
                    break;
                } else {
                    // did not converge: Try lowering the accuracy:
                    tol.val = tol.val * 10;
                    System.out.println("New reduced tolerance: "+tol.val);
                    return solve(eigenvalues, ritz, tol.val, ncvModification);
                }
            }
            if (ido.val != -1 && ido.val != 1)
                throw new IllegalStateException("ido = " + ido.val);
            // could be refactored to handle the other types of mode
            av(workd, ipntr[0] - 1, ipntr[1] - 1);
        }

        log.fine(i + " iterations for " + n);

        if (info.val != 0 && info.val != 1) {
            if (info.val == 3) {
                //  3: No shifts could be applied during a cycle of the
                //     Implicitly restarted Arnoldi iteration. One possibility
                //     is to increase the size of NCV relative to NEV.
                System.err.println("'No shifts could be applied during a cycle!'\n" +
                        "Adding one to the ncv count: "+(ncvModification+1));
                return solve(eigenvalues, ritz, tolerance, ncvModification+1);
            }
            throw new IllegalStateException("info = " + info.val);
        }
        if (info.val == 1){
            // not converged in max num of iterations!
            System.err.println("Maximum number ("+MAX_ITERATIONS+") taken. "+iparam[5]+" converged Ritz values.");
        }
        System.out.println(iparam[5]+" converged 'Ritz' values.");
        double[] dr = new double[nev.val+1];
        double[] di = java.util.Arrays.copyOf(dr,dr.length);
        boolean[] select = new boolean[ncv];
        double[] z = java.util.Arrays.copyOfRange(v, 0, (nev.val+1) * n);

        boolean computeRitzVectors = !this.computeOnlyEigenvalues; // RVEC    LOGICAL  (INPUT)
        // Specifies whether a basis for the invariant subspace corresponding
        // to the converged Ritz value approximations for the eigenproblem
        // A*z = lambda*B*z is computed.
        //            RVEC = .FALSE.     Compute Ritz values only.
        //            RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.

        arpack.dneupd(computeRitzVectors, "P", select, dr, di, z, 1, 0, 0, workd, bmat, n, which, nev, tol.val,
                resid, ncv, v, n, iparam, ipntr, workd, workl, workl.length,
                info);
        if (info.val != 0) {
            System.out.println(Arrays.toString(iparam));
            throw new IllegalStateException("info = " + info.val);
        }

        int computed = iparam[4];
        ArpackGen.log.fine("computed " + computed + " eigenvalues");

        Map<Double, DenseVectorSub> solution = new TreeMap<Double, DenseVectorSub>(
                new Comparator<Double>() {
                    @Override
                    public int compare(Double o1, Double o2) {
                        // highest first
                        return Double.compare(o2, o1);
                    }
                });
        DenseVector eigenvectors;
//        if (iparam[5] > 0) {
//            eigenvectors = new DenseVector(v, false);
//        } else {
        eigenvectors = new DenseVector(v, false);
//        }

        for (i = 0; i < computed; i++) {
            if (new Double(di[i]).compareTo(0.0) == 0) { // ignore complex solutions
                double eigenvalue = dr[i];
                DenseVectorSub eigenvector = new DenseVectorSub(eigenvectors,
                        i * n, n);
                solution.put(eigenvalue, eigenvector);
            }
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
