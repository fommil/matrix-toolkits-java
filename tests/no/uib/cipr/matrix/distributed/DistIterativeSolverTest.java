package no.uib.cipr.matrix.distributed;

import java.util.Arrays;

import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.distributed.CollectiveCommunications;
import no.uib.cipr.matrix.distributed.Communicator;
import no.uib.cipr.matrix.distributed.DistColMatrix;
import no.uib.cipr.matrix.distributed.DistRowMatrix;
import no.uib.cipr.matrix.distributed.DistVector;
import no.uib.cipr.matrix.sparse.BiCGstab;
import no.uib.cipr.matrix.sparse.CG;
import no.uib.cipr.matrix.sparse.DefaultIterationMonitor;
import no.uib.cipr.matrix.sparse.GMRES;
import no.uib.cipr.matrix.sparse.IterationMonitor;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.Utilities;

/**
 * @deprecated the <code>no.uib.cipr.matrix.distributed</code> package has been deprecated because
 * of a number of hard to fix concurrency bugs. It is distributed only for backwards compatibility,
 * but is not recommended. The utility of this package is questionable, as it does not allow
 * distribution of computation between JVMs or across a network. For many people, distributed
 * computing of multiple matrices can be achieved at a user-level through the
 * <a href="http://jppf.org">JPPF Framework</a>.
 * Users who need to deal with few very large matrices may wish to implement their own storage classes
 * and solvers using JPPF, but this will not be supported directly in matrix-toolkits-java.
 */
@Deprecated
public class DistIterativeSolverTest extends TestCase {

    volatile CollectiveCommunications coll;

    volatile DenseMatrix A_unsymm;

    volatile UpperSymmDenseMatrix A_symm;

	volatile DenseVector x, b_unsymm, b_symm;

    /**
     * Partitioning
     */
	volatile int[] localLength;

	volatile double[] output;

    @Override
    protected void setUp() throws Exception {
        int size = Utilities.getInt(1, 8);
        coll = new CollectiveCommunications(size);

        int n = Utilities.getInt(size, 250);
        A_unsymm = new DenseMatrix(n, n);
        A_symm = new UpperSymmDenseMatrix(n);

        Utilities.populate(A_unsymm);
        Utilities.upperPopulate(A_unsymm);

        double shift = 10;

        do {
            Utilities.addDiagonal(A_unsymm, shift);
        } while (Utilities.singular(A_unsymm));

        do {
            Utilities.addDiagonal(A_symm, shift);
        } while (!Utilities.spd(A_symm));

        x = new DenseVector(n);
        b_unsymm = x.copy();
        b_symm = x.copy();

        Utilities.populate(x);
        A_unsymm.mult(x, b_unsymm);
        A_symm.mult(x, b_symm);

        output = new double[n];

        // Set local lengths
        localLength = new int[size];
        Arrays.fill(localLength, n / size);

        // Adjust the last length to ensure the whole vector is covered
        int sum = n;
        for (int l : localLength)
            sum -= l;
        localLength[size - 1] += sum;
    }

    public void testRowGMRES_1() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new GMRESRowDistIterativeSolver(i,
                    Vector.Norm.One));

        compare(t);
    }

    public void testRowGMRES_2() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new GMRESRowDistIterativeSolver(i,
                    Vector.Norm.Two));

        compare(t);
    }

    public void testRowGMRES_inf() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new GMRESRowDistIterativeSolver(i,
                    Vector.Norm.Infinity));

        compare(t);
    }

    public void testRowBiCGstab_1() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new BiCGstabRowDistIterativeSolver(i,
                    Vector.Norm.One));

        compare(t);
    }

    public void testRowBiCGstab_2() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new BiCGstabRowDistIterativeSolver(i,
                    Vector.Norm.Two));

        compare(t);
    }

    public void testRowBiCGstab_inf() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new BiCGstabRowDistIterativeSolver(i,
                    Vector.Norm.Infinity));

        compare(t);
    }

    public void testColumnGMRES_1() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new GMRESColumnDistIterativeSolver(i,
                    Vector.Norm.One));

        compare(t);
    }

    public void testColumnGMRES_2() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new GMRESColumnDistIterativeSolver(i,
                    Vector.Norm.Two));

        compare(t);
    }

    public void testColumnGMRES_inf() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new GMRESColumnDistIterativeSolver(i,
                    Vector.Norm.Infinity));

        compare(t);
    }

    public void testColumnBiCGstab_1() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new BiCGstabColumnDistIterativeSolver(i,
                    Vector.Norm.One));

        compare(t);
    }

    public void testColumnBiCGstab_2() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new BiCGstabColumnDistIterativeSolver(i,
                    Vector.Norm.Two));

        compare(t);
    }

    public void testColumnBiCGstab_inf() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new BiCGstabColumnDistIterativeSolver(i,
                    Vector.Norm.Infinity));

        compare(t);
    }

    public void testRowCG_1() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new CGRowDistIterativeSolver(i, Vector.Norm.One));

        compare(t);
    }

    public void testRowCG_2() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new CGRowDistIterativeSolver(i, Vector.Norm.Two));

        compare(t);
    }

    public void testRowCG_inf() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new CGRowDistIterativeSolver(i,
                    Vector.Norm.Infinity));

        compare(t);
    }

    public void testColumnCG_1() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new CGColumnDistIterativeSolver(i,
                    Vector.Norm.One));

        compare(t);
    }

    public void testColumnCG_2() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new CGColumnDistIterativeSolver(i,
                    Vector.Norm.Two));

        compare(t);
    }

    public void testColumnCG_inf() throws InterruptedException {
        Thread[] t = new Thread[coll.size()];
        for (int i = 0; i < t.length; ++i)
            t[i] = new Thread(new CGColumnDistIterativeSolver(i,
                    Vector.Norm.Infinity));

        compare(t);
    }

    private void compare(Thread[] t) throws InterruptedException {
		ExecutorService pool =
			Executors.newFixedThreadPool(t.length);
		
        for (Thread ti : t)
            pool.execute(ti);
		
		pool.shutdown();
		pool.awaitTermination(20, TimeUnit.SECONDS);

        for (int i = 0; i < x.size(); ++i)
            assertEquals(x.get(i), output[i], 1e-10);
    }

    private abstract class DistIterativeSolver implements Runnable {
        final protected int rank;

        final protected Vector.Norm norm;

        public DistIterativeSolver(int rank, Vector.Norm norm) {
            this.rank = rank;
            this.norm = norm;
        }

        protected abstract Matrix createMatrix(Communicator comm);

        protected abstract int[] getRowOwnerships(Matrix A);

        protected abstract int[] getColumnOwnerships(Matrix A);

        protected abstract void populateMatrix(Matrix A);

        protected abstract double getVectorEntry(int i);

        protected abstract IterativeSolver createSolver(Vector x);

        public void run() {
            Communicator comm = coll.createCommunicator(rank);

            Matrix A = createMatrix(comm);
            populateMatrix(A);

            DenseVector bl = new DenseVector(localLength[rank]);
            DistVector b_dist = new DistVector(x.size(), comm, bl);
            DistVector x_dist = b_dist.copy();

            int[] n = getRowOwnerships(A);
            for (int i = n[rank]; i < n[rank + 1]; ++i)
                b_dist.set(i, getVectorEntry(i));

            IterativeSolver solver = createSolver(b_dist);

            IterationMonitor monitor = new DefaultIterationMonitor(1000, 1e-50,
                    1e-12, 1e+5);
            monitor.setNormType(norm);
            solver.setIterationMonitor(monitor);

            try {
                solver.solve(A, b_dist, x_dist);
            } catch (IterativeSolverNotConvergedException e) {
                // This will just lead to an error later on
            }

            for (int i = n[rank]; i < n[rank + 1]; ++i)
                output[i] = x_dist.get(i);
        }
    }

    private abstract class RowDistIterativeSolver extends DistIterativeSolver {

        public RowDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected Matrix createMatrix(Communicator comm) {
            int n = x.size();

            Matrix Al = new DenseMatrix(localLength[rank], localLength[rank]);
            Matrix Bl = new DenseMatrix(localLength[rank], n);

            return new DistRowMatrix(n, n, comm, Al, Bl);
        }

        @Override
        protected int[] getColumnOwnerships(Matrix A) {
            return ((DistRowMatrix) A).getColumnOwnerships();
        }

        @Override
        protected int[] getRowOwnerships(Matrix A) {
            return ((DistRowMatrix) A).getRowOwnerships();
        }
    }

    private abstract class ColumnDistIterativeSolver extends
            DistIterativeSolver {

        public ColumnDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected Matrix createMatrix(Communicator comm) {
            int n = x.size();

            Matrix Al = new DenseMatrix(localLength[rank], localLength[rank]);
            Matrix Bl = new DenseMatrix(n, localLength[rank]);

            return new DistColMatrix(n, n, comm, Al, Bl);
        }

        @Override
        protected int[] getColumnOwnerships(Matrix A) {
            return ((DistColMatrix) A).getColumnOwnerships();
        }

        @Override
        protected int[] getRowOwnerships(Matrix A) {
            return ((DistColMatrix) A).getRowOwnerships();
        }
    }

    private abstract class SymmRowDistIterativeSolver extends
            RowDistIterativeSolver {

        public SymmRowDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected void populateMatrix(Matrix A) {
            int[] n = getRowOwnerships(A);
            for (int i = n[rank]; i < n[rank + 1]; ++i)
                for (int j = 0; j < A.numColumns(); ++j)
                    A.set(i, j, A_symm.get(i, j));
        }

        @Override
        protected double getVectorEntry(int i) {
            return b_symm.get(i);
        }
    }

    private abstract class UnSymmRowDistIterativeSolver extends
            RowDistIterativeSolver {

        public UnSymmRowDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected void populateMatrix(Matrix A) {
            int[] n = getRowOwnerships(A);
            for (int i = n[rank]; i < n[rank + 1]; ++i)
                for (int j = 0; j < A.numColumns(); ++j)
                    A.set(i, j, A_unsymm.get(i, j));
        }

        @Override
        protected double getVectorEntry(int i) {
            return b_unsymm.get(i);
        }
    }

    private abstract class SymmColumnDistIterativeSolver extends
            ColumnDistIterativeSolver {

        public SymmColumnDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected void populateMatrix(Matrix A) {
            int[] m = getColumnOwnerships(A);
            for (int i = 0; i < A.numRows(); ++i)
                for (int j = m[rank]; j < m[rank + 1]; ++j)
                    A.set(i, j, A_symm.get(i, j));
        }

        @Override
        protected double getVectorEntry(int i) {
            return b_symm.get(i);
        }
    }

    private abstract class UnSymmColumnDistIterativeSolver extends
            ColumnDistIterativeSolver {

        public UnSymmColumnDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected void populateMatrix(Matrix A) {
            int[] m = getColumnOwnerships(A);
            for (int i = 0; i < A.numRows(); ++i)
                for (int j = m[rank]; j < m[rank + 1]; ++j)
                    A.set(i, j, A_unsymm.get(i, j));
        }

        @Override
        protected double getVectorEntry(int i) {
            return b_unsymm.get(i);
        }
    }

    private class GMRESRowDistIterativeSolver extends
            UnSymmRowDistIterativeSolver {

        public GMRESRowDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected IterativeSolver createSolver(Vector x) {
            return new GMRES(x);
        }
    }

    private class BiCGstabRowDistIterativeSolver extends
            UnSymmRowDistIterativeSolver {

        public BiCGstabRowDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected IterativeSolver createSolver(Vector x) {
            return new BiCGstab(x);
        }
    }

    private class GMRESColumnDistIterativeSolver extends
            UnSymmColumnDistIterativeSolver {

        public GMRESColumnDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected IterativeSolver createSolver(Vector x) {
            return new GMRES(x);
        }
    }

    private class BiCGstabColumnDistIterativeSolver extends
            UnSymmColumnDistIterativeSolver {

        public BiCGstabColumnDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected IterativeSolver createSolver(Vector x) {
            return new BiCGstab(x);
        }
    }

    private class CGRowDistIterativeSolver extends SymmRowDistIterativeSolver {

        public CGRowDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected IterativeSolver createSolver(Vector x) {
            return new CG(x);
        }
    }

    private class CGColumnDistIterativeSolver extends
            SymmColumnDistIterativeSolver {

        public CGColumnDistIterativeSolver(int rank, Vector.Norm norm) {
            super(rank, norm);
        }

        @Override
        protected IterativeSolver createSolver(Vector x) {
            return new CG(x);
        }
    }
}
