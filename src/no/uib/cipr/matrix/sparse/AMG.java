/*
 * Copyright (C) 2003-2006 Bjørn-Ove Heimsund
 * 
 * This file is part of MTJ.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package no.uib.cipr.matrix.sparse;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import no.uib.cipr.matrix.DenseLU;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 * Algebraic multigrid preconditioner. Uses the smoothed aggregation method
 * described by Vanek, Mandel, and Brezina (1996).
 */
public class AMG implements Preconditioner {

    /**
     * Relaxations at each level
     */
    private SSOR[] preM, postM;

    /**
     * The number of levels
     */
    private int m;

    /**
     * System matrix at each level, except at the coarsest
     */
    private CompRowMatrix[] A;

    /**
     * LU factorization at the coarsest level
     */
    private DenseLU lu;

    /**
     * Solution, right-hand side, and residual vectors at each level
     */
    private DenseVector[] u, f, r;

    /**
     * Interpolation operators going to a finer mesh
     */
    private CompColMatrix[] I;

    /**
     * Smallest matrix size before terminating the AMG setup phase. Matrices
     * smaller than this will be solved by a direct solver
     */
    private final int min;

    /**
     * Number of times to perform the pre- and post-smoothings
     */
    private final int nu1, nu2;

    /**
     * Determines cycle type. gamma=1 is V, gamma=2 is W
     */
    private final int gamma;

    /**
     * Overrelaxation parameters in the pre- and post-smoothings, and with the
     * possibility of distinct values in the forward and reverse sweeps
     */
    private final double omegaPreF, omegaPreR, omegaPostF, omegaPostR;

    /**
     * Perform a reverse (backwards) smoothing sweep
     */
    private final boolean reverse;

    /**
     * Jacobi damping parameter, between zero and one. If it equals zero, the
     * method reduces to the standard aggregate multigrid method
     */
    private final double omega;

    /**
     * Operating in transpose mode?
     */
    private boolean transpose;

    /**
     * Sets up the algebraic multigrid preconditioner
     * 
     * @param omegaPreF
     *            Overrelaxation parameter in the forward sweep of the
     *            pre-smoothing
     * @param omegaPreR
     *            Overrelaxation parameter in the backwards sweep of the
     *            pre-smoothing
     * @param omegaPostF
     *            Overrelaxation parameter in the forward sweep of the
     *            post-smoothing
     * @param omegaPostR
     *            Overrelaxation parameter in the backwards sweep of the
     *            post-smoothing
     * @param nu1
     *            Number of pre-relaxations to perform
     * @param nu2
     *            Number of post-relaxations to perform
     * @param gamma
     *            Number of times to go to a coarser level
     * @param min
     *            Smallest matrix size before using a direct solver
     * @param omega
     *            Jacobi damping parameter, between zero and one. If it equals
     *            zero, the method reduces to the standard aggregate multigrid
     *            method
     */
    public AMG(double omegaPreF, double omegaPreR, double omegaPostF,
            double omegaPostR, int nu1, int nu2, int gamma, int min,
            double omega) {
        this.omegaPreF = omegaPreF;
        this.omegaPreR = omegaPreR;
        this.omegaPostF = omegaPostF;
        this.omegaPostR = omegaPostR;

        reverse = true;

        this.nu1 = nu1;
        this.nu2 = nu2;
        this.gamma = gamma;
        this.min = min;

        this.omega = omega;
    }

    /**
     * Sets up the algebraic multigrid preconditioner. Uses an SOR method,
     * without the backward sweep in SSOR
     * 
     * @param omegaPre
     *            Overrelaxation parameter in the pre-smoothing
     * @param omegaPost
     *            Overrelaxation parameter in the post-smoothing
     * @param nu1
     *            Number of pre-relaxations to perform
     * @param nu2
     *            Number of post-relaxations to perform
     * @param gamma
     *            Number of times to go to a coarser level
     * @param min
     *            Smallest matrix size before using a direct solver
     * @param omega
     *            Jacobi damping parameter, between zero and one. If it equals
     *            zero, the method reduces to the standard aggregate multigrid
     *            method
     */
    public AMG(double omegaPre, double omegaPost, int nu1, int nu2, int gamma,
            int min, double omega) {
        this.omegaPreF = omegaPre;
        this.omegaPreR = omegaPre;
        this.omegaPostF = omegaPost;
        this.omegaPostR = omegaPost;

        reverse = false;

        this.nu1 = nu1;
        this.nu2 = nu2;
        this.gamma = gamma;
        this.min = min;

        this.omega = omega;
    }

    /**
     * Sets up the algebraic multigrid preconditioner using some default
     * parameters. In the presmoothing, <code>omegaF=1</code> and
     * <code>omegaR=1.85</code>, while in the postsmoothing,
     * <code>omegaF=1.85</code> and <code>omegaR=1</code>. Sets
     * <code>nu1=nu2=gamma=1</code>, has a smallest matrix size of 40, and
     * sets <code>omega=2/3</code>.
     */
    public AMG() {
        this(1, 1.85, 1.85, 1, 1, 1, 1, 40, 2. / 3);
    }

    public Vector apply(Vector b, Vector x) {
        u[0].set(x);
        f[0].set(b);

        transpose = false;
        cycle(0);

        return x.set(u[0]);
    }

    public Vector transApply(Vector b, Vector x) {
        u[0].set(x);
        f[0].set(b);

        transpose = true;
        cycle(0);

        return x.set(u[0]);
    }

    public void setMatrix(Matrix A) {
        List<CompRowMatrix> Al = new LinkedList<CompRowMatrix>();
        List<CompColMatrix> Il = new LinkedList<CompColMatrix>();

        Al.add(new CompRowMatrix(A));

        for (int k = 0; Al.get(k).numRows() > min; ++k) {

            CompRowMatrix Af = Al.get(k);

            double eps = 0.08 * Math.pow(0.5, k);

            // Create the aggregates
            Aggregator aggregator = new Aggregator(Af, eps);

            // If no aggregates were created, no interpolation operator will be
            // created, and the setup phase stops
            if (aggregator.getAggregates().size() == 0)
                break;

            // Create an interpolation operator using smoothing. This also
            // creates the Galerkin operator
            Interpolator sa = new Interpolator(aggregator, Af, omega);

            Al.add(sa.getGalerkinOperator());
            Il.add(sa.getInterpolationOperator());
        }

        // Copy to array storage
        m = Al.size();
        if (m == 0)
            throw new RuntimeException("Matrix too small for AMG");

        I = new CompColMatrix[m - 1];
        this.A = new CompRowMatrix[m - 1];

        Il.toArray(I);
        for (int i = 0; i < Al.size() - 1; ++i)
            this.A[i] = Al.get(i);

        // Create a LU decomposition of the smallest Galerkin matrix
        DenseMatrix Ac = new DenseMatrix(Al.get(Al.size() - 1));
        lu = new DenseLU(Ac.numRows(), Ac.numColumns());
        lu.factor(Ac);

        // Allocate vectors at each level
        u = new DenseVector[m];
        f = new DenseVector[m];
        r = new DenseVector[m];
        for (int k = 0; k < m; ++k) {
            int n = Al.get(k).numRows();
            u[k] = new DenseVector(n);
            f[k] = new DenseVector(n);
            r[k] = new DenseVector(n);
        }

        // Set up the SSOR relaxation schemes
        preM = new SSOR[m - 1];
        postM = new SSOR[m - 1];
        for (int k = 0; k < m - 1; ++k) {
            CompRowMatrix Ak = this.A[k];
            preM[k] = new SSOR(Ak, reverse, omegaPreF, omegaPreR);
            postM[k] = new SSOR(Ak, reverse, omegaPostF, omegaPostR);
            preM[k].setMatrix(Ak);
            postM[k].setMatrix(Ak);
        }
    }

    /**
     * Performs a multigrid cycle
     * 
     * @param k
     *            Level to cycle at. Start by calling <code>cycle(0)</code>
     */
    private void cycle(int k) {
        if (k == m - 1)
            directSolve();
        else {

            // Presmoothings
            preRelax(k);

            u[k + 1].zero();

            // Compute the residual
            A[k].multAdd(-1, u[k], r[k].set(f[k]));

            // Restrict to the next coarser level
            I[k].transMult(r[k], f[k + 1]);

            // Recurse to next level
            for (int i = 0; i < gamma; ++i)
                cycle(k + 1);

            // Add residual correction by prolongation
            I[k].multAdd(u[k + 1], u[k]);

            // Postsmoothings
            postRelax(k);
        }
    }

    /**
     * Solves directly at the coarsest level
     */
    private void directSolve() {
        int k = m - 1;
        u[k].set(f[k]);
        DenseMatrix U = new DenseMatrix(u[k], false);

        if (transpose)
            lu.transSolve(U);
        else
            lu.solve(U);
    }

    /**
     * Applies the relaxation scheme at the given level
     * 
     * @param k
     *            Multigrid level
     */
    private void preRelax(int k) {
        for (int i = 0; i < nu1; ++i)
            if (transpose)
                preM[k].transApply(f[k], u[k]);
            else
                preM[k].apply(f[k], u[k]);
    }

    /**
     * Applies the relaxation scheme at the given level
     * 
     * @param k
     *            Multigrid level
     */
    private void postRelax(int k) {
        for (int i = 0; i < nu2; ++i)
            if (transpose)
                postM[k].transApply(f[k], u[k]);
            else
                postM[k].apply(f[k], u[k]);
    }

    /**
     * Creates aggregates. These are disjoint sets, each of which represents one
     * node at a coarser mesh by aggregating together a set of fine nodes
     */
    private static class Aggregator {

        /**
         * The aggregates
         */
        private List<Set<Integer>> C;

        /**
         * Diagonal indices into the sparse matrix
         */
        private int[] diagind;

        /**
         * The strongly coupled node neighborhood of a given node
         */
        private List<Set<Integer>> N;

        /**
         * Creates the aggregates
         * 
         * @param A
         *            Sparse matrix
         * @param eps
         *            Tolerance for selecting the strongly coupled node
         *            neighborhoods. Between zero and one.
         */
        public Aggregator(CompRowMatrix A, double eps) {

            diagind = findDiagonalIndices(A);
            N = findNodeNeighborhood(A, diagind, eps);

            /*
             * Initialization. Remove isolated nodes from the aggregates
             */

            boolean[] R = createInitialR(A);

            /*
             * Startup aggregation. Use disjoint strongly coupled neighborhoods
             * as the initial aggregate approximation
             */

            C = createInitialAggregates(N, R);

            /*
             * Enlargment of the aggregates. Add nodes to each aggregate based
             * on how strongly connected the nodes are to a given aggregate
             */

            C = enlargeAggregates(C, N, R);

            /*
             * Handling of the remenants. Put all remaining unallocated nodes
             * into new aggregates defined by the intersection of N and R
             */

            C = createFinalAggregates(C, N, R);
        }

        /**
         * Gets the aggregates
         */
        public List<Set<Integer>> getAggregates() {
            return C;
        }

        /**
         * Returns the matrix diagonal indices. This is a by-product of the
         * aggregation
         */
        public int[] getDiagonalIndices() {
            return diagind;
        }

        /**
         * Returns the strongly coupled node neighborhoods of a given node. This
         * is a by-product of the aggregation
         */
        public List<Set<Integer>> getNodeNeighborhoods() {
            return N;
        }

        /**
         * Finds the diagonal indices of the matrix
         */
        private int[] findDiagonalIndices(CompRowMatrix A) {
            int[] rowptr = A.getRowPointers();
            int[] colind = A.getColumnIndices();

            int[] diagind = new int[A.numRows()];

            for (int i = 0; i < A.numRows(); ++i) {
                diagind[i] = no.uib.cipr.matrix.sparse.Arrays.binarySearch(
                        colind, i, rowptr[i], rowptr[i + 1]);
                if (diagind[i] < 0)
                    throw new RuntimeException(
                            "Matrix is missing a diagonal entry on row "
                                    + (i + 1));
            }

            return diagind;
        }

        /**
         * Finds the strongly coupled node neighborhoods
         */
        private List<Set<Integer>> findNodeNeighborhood(CompRowMatrix A,
                int[] diagind, double eps) {

            N = new ArrayList<Set<Integer>>(A.numRows());

            int[] rowptr = A.getRowPointers();
            int[] colind = A.getColumnIndices();
            double[] data = A.getData();

            for (int i = 0; i < A.numRows(); ++i) {
                Set<Integer> Ni = new HashSet<Integer>();

                double aii = data[diagind[i]];
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {
                    double aij = data[j];
                    double ajj = data[diagind[colind[j]]];

                    if (Math.abs(aij) >= eps * Math.sqrt(aii * ajj))
                        Ni.add(colind[j]);
                }

                N.add(Ni);
            }

            return N;
        }

        /**
         * Creates the initial R-set by including only the connected nodes
         */
        private boolean[] createInitialR(CompRowMatrix A) {
            boolean[] R = new boolean[A.numRows()];

            int[] rowptr = A.getRowPointers();
            int[] colind = A.getColumnIndices();
            double[] data = A.getData();

            for (int i = 0; i < A.numRows(); ++i) {
                boolean hasOffDiagonal = false;

                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
                    if (colind[j] != i && data[j] != 0) {
                        hasOffDiagonal = true;
                        break;
                    }

                R[i] = hasOffDiagonal;
            }

            return R;
        }

        /**
         * Creates the initial aggregates
         */
        private List<Set<Integer>> createInitialAggregates(
                List<Set<Integer>> N, boolean[] R) {
            C = new ArrayList<Set<Integer>>();

            for (int i = 0; i < R.length; ++i) {

                // Skip non-free nodes
                if (!R[i])
                    continue;

                // See if all nodes in the current N-set are free
                boolean free = true;
                for (int j : N.get(i))
                    free &= R[j];

                // Create an aggregate out of N[i]
                if (free) {
                    C.add(new HashSet<Integer>(N.get(i)));
                    for (int j : N.get(i))
                        R[j] = false;
                }

            }

            return C;
        }

        /**
         * Enlarges the aggregates
         */
        private List<Set<Integer>> enlargeAggregates(List<Set<Integer>> C,
                List<Set<Integer>> N, boolean[] R) {

            // Contains the aggregates each node is coupled to
            List<List<Integer>> belong = new ArrayList<List<Integer>>(R.length);
            for (int i = 0; i < R.length; ++i)
                belong.add(new ArrayList<Integer>());

            // Find which aggregate each node is coupled to. This is used for
            // the intersection between Ni and Ck
            for (int k = 0; k < C.size(); ++k)
                for (int j : C.get(k))
                    belong.get(j).add(k);

            // Number of nodes in the intersection between each C and Ni
            int[] intersect = new int[C.size()];

            for (int i = 0; i < R.length; ++i) {

                // Skip non-free nodes
                if (!R[i])
                    continue;

                // Find the number of nodes intersecting Ni and every C, and
                // keep a track on the largest overlap
                Arrays.fill(intersect, 0);
                int largest = 0, maxValue = 0;
                for (int j : N.get(i))

                    // The k-index is to an aggregate coupled to node j
                    for (int k : belong.get(j)) {
                        intersect[k]++;
                        if (intersect[k] > maxValue) {
                            largest = k;
                            maxValue = intersect[largest];
                        }
                    }

                // Add the node to the proper C-set, and mark it as used
                // Also, check if the node actually does couple to a set
                if (maxValue > 0) {
                    R[i] = false;
                    C.get(largest).add(i);
                }
            }

            return C;
        }

        /**
         * Creates final aggregates from the remaining unallocated nodes
         */
        private List<Set<Integer>> createFinalAggregates(List<Set<Integer>> C,
                List<Set<Integer>> N, boolean[] R) {

            for (int i = 0; i < R.length; ++i) {

                // Skip non-free nodes
                if (!R[i])
                    continue;

                // Create new aggregate from the nodes in N[i] which are free
                Set<Integer> Cn = new HashSet<Integer>();
                for (int j : N.get(i))
                    if (R[j]) {
                        R[j] = false;
                        Cn.add(j);
                    }

                if (!Cn.isEmpty())
                    C.add(Cn);
            }

            return C;
        }
    }

    /**
     * Creates interpolation (prolongation) operators using based on the
     * aggregates. Can optionally smooth the aggregates
     */
    private static class Interpolator {

        /**
         * The Galerkin coarse-space operator
         */
        private CompRowMatrix Ac;

        /**
         * The interpolation (prolongation) matrix
         */
        private CompColMatrix I;

        /**
         * Creates the interpolation (prolongation) and Galerkin operators
         * 
         * @param aggregator
         *            Aggregates
         * @param A
         *            Matrix
         * @param omega
         *            Jacobi damping parameter between zero and one. If zero, no
         *            smoothing is performed, and a faster algorithm for forming
         *            the Galerkin operator will be used.
         */
        public Interpolator(Aggregator aggregator, CompRowMatrix A, double omega) {
            List<Set<Integer>> C = aggregator.getAggregates();
            List<Set<Integer>> N = aggregator.getNodeNeighborhoods();
            int[] diagind = aggregator.getDiagonalIndices();

            // Create the tentative prolongation, in compressed form
            int[] pt = createTentativeProlongation(C, A.numRows());

            /*
             * Apply Jacobi smoothing to the prolongator
             */

            if (omega != 0) {

                // Smooth the operator by a damped Jacobi method
                List<Map<Integer, Double>> P = createSmoothedProlongation(C, N,
                        A, diagind, omega, pt);

                // Form a compressed column storage for the operator
                I = createInterpolationMatrix(P, A.numRows());

                // Create the Galerkin operator using a slow method
                Ac = createGalerkinSlow(I, A);
            }

            /*
             * Use the aggregates as-is
             */

            else {

                // Create the Galerkin operator using a fast method
                Ac = createGalerkinFast(A, pt, C.size());

                // Form an explicit interpolation operator
                I = createInterpolationMatrix(pt, C.size());
            }
        }

        /**
         * Creates the tentative prolongation operator. Since the columns are
         * all disjoint, and its entries are binary, it is possible to store it
         * in a single array. Its length equals the number of fine nodes, and
         * the entries are the indices to the corresponding aggregate (C-set).
         */
        private int[] createTentativeProlongation(List<Set<Integer>> C, int n) {
            int[] pt = new int[n];
            Arrays.fill(pt, -1);

            for (int i = 0; i < C.size(); ++i)
                for (int j : C.get(i))
                    pt[j] = i;

            return pt;
        }

        /**
         * Creates the Galerkin operator using the assumption of disjoint
         * (non-smoothed) aggregates
         */
        private CompRowMatrix createGalerkinFast(CompRowMatrix A, int[] pt,
                int c) {
            int n = pt.length;

            FlexCompRowMatrix Ac = new FlexCompRowMatrix(c, c);

            int[] rowptr = A.getRowPointers();
            int[] colind = A.getColumnIndices();
            double[] data = A.getData();

            for (int i = 0; i < n; ++i)
                if (pt[i] != -1)
                    for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
                        if (pt[colind[j]] != -1)
                            Ac.add(pt[i], pt[colind[j]], data[j]);

            return new CompRowMatrix(Ac);
        }

        /**
         * Creates the interpolation (prolongation) matrix based on the smoothed
         * aggregates
         */
        private CompColMatrix createInterpolationMatrix(
                List<Map<Integer, Double>> P, int n) {

            // Determine the sparsity pattern of I
            int c = P.size();
            int[][] nz = new int[c][];
            for (int j = 0; j < c; ++j) {

                Map<Integer, Double> Pj = P.get(j);
                nz[j] = new int[Pj.size()];

                int l = 0;
                for (int k : Pj.keySet())
                    nz[j][l++] = k;
            }

            I = new CompColMatrix(n, c, nz);

            // Populate it with numerical entries
            for (int j = 0; j < c; ++j) {

                Map<Integer, Double> Pj = P.get(j);

                for (Map.Entry<Integer, Double> e : Pj.entrySet())
                    I.set(e.getKey(), j, e.getValue());
            }

            return I;
        }

        /**
         * Creates the interpolation (prolongation) matrix based on the
         * non-smoothed aggregates
         */
        private CompColMatrix createInterpolationMatrix(int[] pt, int c) {
            FlexCompColMatrix If = new FlexCompColMatrix(pt.length, c);

            for (int i = 0; i < pt.length; ++i)
                if (pt[i] != -1)
                    If.set(i, pt[i], 1);

            return new CompColMatrix(If);
        }

        /**
         * Gets the interpolation (prolongation) operator
         */
        public CompColMatrix getInterpolationOperator() {
            return I;
        }

        /**
         * Creates the smoothes interpolation (prolongation) operator by a
         * single sweep of the damped Jacobi method
         */
        private List<Map<Integer, Double>> createSmoothedProlongation(
                List<Set<Integer>> C, List<Set<Integer>> N, CompRowMatrix A,
                int[] diagind, double omega, int[] pt) {

            int n = A.numRows(), c = C.size();

            // Allocate the interpolation (prolongation) operator
            // It is stored by columns, so the maps take row-indices as keys
            List<Map<Integer, Double>> P = new ArrayList<Map<Integer, Double>>(
                    c);
            for (int i = 0; i < c; ++i)
                P.add(new HashMap<Integer, Double>());

            int[] rowptr = A.getRowPointers();
            int[] colind = A.getColumnIndices();
            double[] data = A.getData();

            double[] dot = new double[c];

            // Apply the damped Jacobi smoother
            for (int i = 0; i < n; ++i) {

                if (pt[i] == -1)
                    continue;

                Arrays.fill(dot, 0);
                Set<Integer> Ni = N.get(i);

                // Calculate A*Pt, except for the diagonal
                double weakAij = 0;
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j) {

                    if (pt[colind[j]] == -1)
                        continue;

                    double aij = data[j];

                    // Off-diagonal, include only strong couplings, and add the
                    // weak couplings to the diagonal
                    if (aij != 0 && !Ni.contains(colind[j])) {
                        weakAij += aij;
                        continue;
                    }

                    dot[pt[colind[j]]] += aij;
                }

                // Subtract the weak couplings from the diagonal part of A*Pt
                dot[pt[i]] -= weakAij;

                // Scale by omega and the inverse of the diagonal (damping)
                double scale = -omega / data[diagind[i]];
                for (int j = 0; j < dot.length; ++j)
                    dot[j] *= scale;

                // Set to (I-omega*D^{-1}*A)*Pt
                dot[pt[i]]++;

                // This has formed a whole row of P=(I-omega*D^{-1}*A)*Pt
                // Store the non-zeros into the sparse structure
                for (int j = 0; j < dot.length; ++j)
                    if (dot[j] != 0)
                        P.get(j).put(i, dot[j]);
            }

            return P;
        }

        /**
         * Creates the entries of the Galerkin operator
         * <code>Ac = I<sup>T</sup> A I</code>. This is a very
         * time-consuming operation
         */
        private CompRowMatrix createGalerkinSlow(CompColMatrix I,
                CompRowMatrix A) {
            int n = I.numRows(), c = I.numColumns();
            FlexCompRowMatrix Ac = new FlexCompRowMatrix(c, c);

            double[] aiCol = new double[n];
            double[] iCol = new double[n];
            DenseVector aiV = new DenseVector(aiCol, false);
            DenseVector iV = new DenseVector(iCol, false);
            double[] itaiCol = new double[c];
            DenseVector itaiV = new DenseVector(itaiCol, false);

            int[] colptr = I.getColumnPointers();
            int[] rowind = I.getRowIndices();
            double[] Idata = I.getData();

            for (int k = 0; k < c; ++k) {

                // Expand column 'k' of I to dense storage
                iV.zero();
                for (int i = colptr[k]; i < colptr[k + 1]; ++i)
                    iCol[rowind[i]] = Idata[i];

                // Form column 'k' of A*I
                A.mult(iV, aiV);

                // Form column 'k' of I'*A*I
                I.transMult(aiV, itaiV);

                // Store non-zeros into Ac
                for (int i = 0; i < c; ++i)
                    if (itaiCol[i] != 0)
                        Ac.set(i, k, itaiCol[i]);
            }

            return new CompRowMatrix(Ac);
        }

        /**
         * Gets the Galerkin operator
         */
        public CompRowMatrix getGalerkinOperator() {
            return Ac;
        }

    }
}
