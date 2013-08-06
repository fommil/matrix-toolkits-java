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

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 * Iterative linear solver. Solves <code>Ax=b</code> for <code>x</code>,
 * and it supports preconditioning and convergence monitoring.
 */
public interface IterativeSolver {

    /**
     * Solves the given problem, writing result into the vector.
     * 
     * @param A
     *            Matrix of the problem
     * @param b
     *            Right hand side
     * @param x
     *            Solution is stored here. Also used as initial guess
     * @return The solution vector x
     */
    Vector solve(Matrix A, Vector b, Vector x)
            throws IterativeSolverNotConvergedException;

    /**
     * Sets preconditioner
     * 
     * @param M
     *            Preconditioner to use
     */
    void setPreconditioner(Preconditioner M);

    /**
     * Gets preconditioner
     * 
     * @return Current preconditioner
     */
    Preconditioner getPreconditioner();

    /**
     * Sets iteration monitor
     * 
     * @param iter
     *            Iteration monitor
     */
    void setIterationMonitor(IterationMonitor iter);

    /**
     * Gets the iteration monitor
     * 
     * @return Current iteration monitor
     */
    IterationMonitor getIterationMonitor();

}
