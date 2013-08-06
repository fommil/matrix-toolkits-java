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
 * Partial implementation of an iterative solver
 */
public abstract class AbstractIterativeSolver implements IterativeSolver {

    /**
     * Preconditioner to use
     */
    protected Preconditioner M;

    /**
     * Iteration monitor
     */
    protected IterationMonitor iter;

    /**
     * Constructor for AbstractIterativeSolver. Does not use preconditioning,
     * and uses the default linear iteration object.
     */
    public AbstractIterativeSolver() {
        M = new IdentityPreconditioner();
        iter = new DefaultIterationMonitor();
    }

    public void setPreconditioner(Preconditioner M) {
        this.M = M;
    }

    public Preconditioner getPreconditioner() {
        return M;
    }

    public IterationMonitor getIterationMonitor() {
        return iter;
    }

    public void setIterationMonitor(IterationMonitor iter) {
        this.iter = iter;
    }

    /**
     * Checks sizes of input data for {@link #solve(Matrix, Vector, Vector)}.
     * Throws an exception if the sizes does not match.
     */
    protected void checkSizes(Matrix A, Vector b, Vector x) {
        if (!A.isSquare())
            throw new IllegalArgumentException("!A.isSquare()");
        if (b.size() != A.numRows())
            throw new IllegalArgumentException("b.size() != A.numRows()");
        if (b.size() != x.size())
            throw new IllegalArgumentException("b.size() != x.size()");
    }

    /**
     * Identity preconditioner which does nothing
     */
    private static class IdentityPreconditioner implements Preconditioner {

        public Vector apply(Vector b, Vector x) {
            return x.set(b);
        }

        public Vector transApply(Vector b, Vector x) {
            return x.set(b);
        }

        public void setMatrix(Matrix A) {
            // nothing to do
        }

    }

}
