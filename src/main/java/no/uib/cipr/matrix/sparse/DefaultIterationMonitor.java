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
 * Default iteration monitor. This tester checks declares convergence if the
 * absolute value of the residual norm is sufficiently small, or if the relative
 * decrease is small. Divergence is decleared if too many iterations are spent,
 * or the residual has grown too much. NaNs will also cause divergence to be
 * flagged.
 */
public class DefaultIterationMonitor extends AbstractIterationMonitor {

    /**
     * Initial residual
     */
    double initR;

    /**
     * Relative tolerance
     */
    double rtol;

    /**
     * Absolute tolerance
     */
    double atol;

    /**
     * Divergence tolerance
     */
    double dtol;

    /**
     * Maximum number of iterations
     */
    int maxIter;

    /**
     * Constructor for DefaultIterationMonitor
     * 
     * @param maxIter
     *            Maximum number of iterations
     * @param rtol
     *            Relative convergence tolerance (to initial residual)
     * @param atol
     *            Absolute convergence tolerance
     * @param dtol
     *            Relative divergence tolerance (to initial residual)
     */
    public DefaultIterationMonitor(int maxIter, double rtol, double atol,
            double dtol) {
        this.maxIter = maxIter;
        this.rtol = rtol;
        this.atol = atol;
        this.dtol = dtol;
    }

    /**
     * Constructor for DefaultIterationMonitor. Default is 100000 iterations at
     * most, relative tolerance of 1e-5, absolute tolerance of 1e-50 and a
     * divergence tolerance of 1e+5.
     */
    public DefaultIterationMonitor() {
        this.maxIter = 100000;
        this.rtol = 1e-5;
        this.atol = 1e-50;
        this.dtol = 1e+5;
    }

    /**
     * Sets maximum number of iterations to permit
     * 
     * @param maxIter
     *            Maximum number of iterations
     */
    public void setMaxIterations(int maxIter) {
        this.maxIter = maxIter;
    }

    /**
     * Sets the relative tolerance
     * 
     * @param rtol
     *            Relative convergence tolerance (to initial residual)
     */
    public void setRelativeTolerance(double rtol) {
        this.rtol = rtol;
    }

    /**
     * Sets the absolute tolerance
     * 
     * @param atol
     *            Absolute convergence tolerance
     */
    public void setAbsoluteTolerance(double atol) {
        this.atol = atol;
    }

    /**
     * Sets the divergence tolerance
     * 
     * @param dtol
     *            Relative divergence tolerance (to initial residual)
     */
    public void setDivergenceTolerance(double dtol) {
        this.dtol = dtol;
    }

    @Override
    protected boolean convergedI(double r)
            throws IterativeSolverNotConvergedException {
        // Store initial residual
        if (isFirst())
            initR = r;

        // Check for convergence
        if (r < Math.max(rtol * initR, atol))
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
    protected boolean convergedI(double r, Vector x)
            throws IterativeSolverNotConvergedException {
        return convergedI(r);
    }

}
