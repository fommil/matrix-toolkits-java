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
 * Preconditioner interface. Before a preconditioner is used,
 * <code>setMatrix</code> must be called
 */
public interface Preconditioner {

    /**
     * Solves the approximate problem with the given right hand side. Result is
     * stored in given solution vector
     * 
     * @param b
     *            Right hand side of problem
     * @param x
     *            Result is stored here
     * @return x
     */
    Vector apply(Vector b, Vector x);

    /**
     * Solves the approximate transpose problem with the given right hand side.
     * Result is stored in given solution vector
     * 
     * @param b
     *            Right hand side of problem
     * @param x
     *            Result is stored here
     * @return x
     */
    Vector transApply(Vector b, Vector x);

    /**
     * Sets the operator matrix for the preconditioner. This method must be
     * called before a preconditioner is used by an iterative solver
     * 
     * @param A
     *            Matrix to setup the preconditioner for. Not modified
     */
    void setMatrix(Matrix A);

}
