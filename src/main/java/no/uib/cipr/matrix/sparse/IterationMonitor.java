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

import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;

/**
 * Monitors the iterative solution process for convergence and divergence. Can
 * also report the current progress.
 */
public interface IterationMonitor {

    /**
     * Resets the iteration
     */
    void setFirst();

    /**
     * Returns true for the first iteration
     */
    boolean isFirst();

    /**
     * Increases iteration counter
     */
    void next();

    /**
     * Number of iterations performed
     */
    int iterations();

    /**
     * Returns current residual
     */
    double residual();

    /**
     * Checks for convergence
     * 
     * @param r
     *            Residual-vector
     * @param x
     *            State-vector
     * @return True if converged
     */
    boolean converged(Vector r, Vector x)
            throws IterativeSolverNotConvergedException;

    /**
     * Checks for convergence
     * 
     * @param r
     *            Residual-norm
     * @param x
     *            State-vector
     * @return True if converged
     */
    boolean converged(double r, Vector x)
            throws IterativeSolverNotConvergedException;

    /**
     * Checks for convergence
     * 
     * @param r
     *            Residual-norm
     * @return True if converged
     */
    boolean converged(double r) throws IterativeSolverNotConvergedException;

    /**
     * Checks for convergence
     * 
     * @param r
     *            Residual-vector
     * @return True if converged
     */
    boolean converged(Vector r) throws IterativeSolverNotConvergedException;

    /**
     * Sets new iteration reporter
     */
    void setIterationReporter(IterationReporter monitor);

    /**
     * Returns current iteration reporter
     */
    IterationReporter getIterationReporter();

    /**
     * Sets the vector-norm to calculate with
     */
    void setNormType(Norm normType);

    /**
     * Returns the vector-norm in use
     */
    Norm getNormType();

}
