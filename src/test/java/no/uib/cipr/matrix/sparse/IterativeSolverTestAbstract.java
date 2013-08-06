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

import no.uib.cipr.matrix.DenseLU;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.sparse.Preconditioner;
import no.uib.cipr.matrix.Utilities;
import junit.framework.TestCase;

/**
 * Test of the iterative solvers and preconditioners
 */
public abstract class IterativeSolverTestAbstract extends TestCase {

    /**
     * Number of times to repeat tests
     */
    private int repeat = 5;

    /**
     * Sizes of the system matrix
     */
    protected int max = 50, bmax = 10;

    /**
     * Numerical tolerance
     */
    protected double tol = 1e-4;

    /**
     * Diagonal shift for singularity handling
     */
    protected double shift = 100;

    /**
     * Square system matrix
     */
    protected Matrix A;

    /**
     * Right hand side, right hand for transpose system, and the solution vector
     * in both cases
     */
    protected Vector b, bt, x;

    /**
     * Stores the data of x
     */
    protected double[] xd;

    /**
     * Iterative solver to use
     */
    protected IterativeSolver solver;

    /**
     * Preconditioner to use
     */
    protected Preconditioner M;

    /**
     * Constructor for IterativeSolverTestAbstract
     */
    public IterativeSolverTestAbstract(String arg0) {
        super(arg0);
    }

    @Override
    protected void setUp() throws Exception {
        createMatrix();

        int n = A.numRows();
        x = Matrices.random(n);
        b = Matrices.random(n);
        bt = Matrices.random(n);

        // Create solver and preconditioner
        createSolver();
        if (M == null)
            M = solver.getPreconditioner();
        M.setMatrix(A);

        // Compute the correct right hand sides
        b = A.mult(x, b);
        bt = A.transMult(x, bt);

        // Store x for later. It is overwritten
        xd = Matrices.getArray(x);

        // Randomize the inital solution vector
        Matrices.random(x);
    }

    protected abstract void createSolver() throws Exception;

    protected void createMatrix() throws Exception {
        // Create an arbitrary matrix
        int n = Utilities.getInt(1, max);
        int b = Utilities.getInt(Math.min(bmax, n));
        A = new FlexCompRowMatrix(n, n);
        Utilities.rowPopulate(A, b);

        // Make it non-singular
        addDiagonal(A, shift);
        DenseLU lu = DenseLU.factorize(A);
        while (lu.isSingular()) {
            addDiagonal(A, shift);
            lu = DenseLU.factorize(A);
        }
    }

    protected void addDiagonal(Matrix A, double shift) {
        int n = A.numRows(), m = A.numColumns();
        for (int i = 0; i < Math.min(n, m); ++i)
            A.add(i, i, shift);
    }

    @Override
    protected void tearDown() throws Exception {
        A = null;
        b = bt = x = null;
        xd = null;
        solver = null;
    }

    public void testSolve() {
        try {
            solver.solve(A, b, x);
            assertSolved();
        } catch (IterativeSolverNotConvergedException e) {
            fail("Solver did not converge: " + e.getReason() + ". Residual="
                    + e.getResidual());
        }
    }

    public void testRepeatSolve() {
        try {
            for (int i = 0; i < repeat; ++i) {
                solver.solve(A, b, x);
                assertSolved();
                x = Matrices.random(A.numRows());
            }
        } catch (IterativeSolverNotConvergedException e) {
            fail("Solver did not converge: " + e.getReason() + ". Residual="
                    + e.getResidual());
        }
    }

    protected void assertSolved() {
        for (int i = 0; i < xd.length; ++i)
            assertEquals(xd[i], x.get(i), tol);
    }

}
