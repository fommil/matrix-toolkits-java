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

import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.Vector;

/**
 * Iteration monitor based on matrix norms. Extends the default linear iteration
 * object to compare with the norm of the system matrix and the right hand side.
 * Can often be a better convergence criteria than the default, but requires the
 * computation of the matrix norm.
 */
public class MatrixIterationMonitor extends DefaultIterationMonitor {

    /**
     * Norm of the system matrix
     */
    private double normA;

    /**
     * Norm of the right hand side
     */
    private double normb;

    /**
     * Constructor for MatrixIterationMonitor
     * 
     * @param normA
     *            Norm of the matrix A
     * @param normb
     *            Norm of the vector b
     * @param maxIter
     *            Maximum number of iterations
     * @param rtol
     *            Relative convergence tolerance (to initial residual)
     * @param atol
     *            Absolute convergence tolerance
     * @param dtol
     *            Relative divergence tolerance (to initial residual)
     */
    public MatrixIterationMonitor(double normA, double normb, int maxIter,
            double rtol, double atol, double dtol) {
        this.normA = normA;
        this.normb = normb;
        this.maxIter = maxIter;
        this.rtol = rtol;
        this.atol = atol;
        this.dtol = dtol;
    }

    /**
     * Constructor for MatrixIterationMonitor. Default is 100000 iterations at
     * most, relative tolerance of 1e-5, absolute tolerance of 1e-50 and a
     * divergence tolerance of 1e+5.
     */
    public MatrixIterationMonitor(double normA, double normb) {
        this.normA = normA;
        this.normb = normb;
    }

    /**
     * Sets the norm of the system matrix
     * 
     * @param normA
     *            Norm of the matrix A
     */
    public void setMatrixNorm(double normA) {
        this.normA = normA;
    }

    /**
     * Sets the norm of the right hand side vector
     * 
     * @param normb
     *            Norm of the vector b
     */
    public void setVectorNorm(double normb) {
        this.normb = normb;
    }

    @Override
    protected boolean convergedI(double r, Vector x)
            throws IterativeSolverNotConvergedException {
        // Store initial residual
        if (isFirst())
            initR = r;

        // Check for convergence
        if (r < Math.max(rtol * (normA * x.norm(normType) + normb), atol))
            return true;

        // Check for divergence
        if (r > dtol * initR)
            throw new IterativeSolverNotConvergedException(
                    NotConvergedException.Reason.Divergence, this);
        if (iter >= maxIter)
            throw new IterativeSolverNotConvergedException(
                    NotConvergedException.Reason.Iterations, this);
        if (Double.isNaN(r))
            throw new IterativeSolverNotConvergedException(
                    NotConvergedException.Reason.Divergence, this);

        // Neither convergence nor divergence
        return false;
    }

    @Override
    protected boolean convergedI(double r) {
        throw new UnsupportedOperationException();
    }

}
