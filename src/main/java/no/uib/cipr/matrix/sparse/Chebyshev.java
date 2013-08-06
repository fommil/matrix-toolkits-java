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
/*
 * Derived from public domain software at http://www.netlib.org/templates
 */

package no.uib.cipr.matrix.sparse;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 * Chebyshev solver. Solves the symmetric positive definite linear system
 * <code>Ax = b</code> using the Preconditioned Chebyshev Method. Chebyshev
 * requires an acurate estimate on the bounds of the spectrum of the matrix.
 * 
 * @author Templates
 */
public class Chebyshev extends AbstractIterativeSolver {

    /**
     * Estimates for the eigenvalue of the matrix
     */
    private double eigmin, eigmax;

    /**
     * Vectors for use in the iterative solution process
     */
    private Vector p, z, r, q;

    /**
     * Constructor for Chebyshev. Uses the given vector as template for creating
     * scratch vectors. Typically, the solution or the right hand side vector
     * can be passed, and the template is not modified. Eigenvalue estimates
     * must also be provided
     * 
     * @param template
     *            Vector to use as template for the work vectors needed in the
     *            solution process
     * @param eigmin
     *            Smallest eigenvalue. Must be positive
     * @param eigmax
     *            Largest eigenvalue. Must be positive
     */
    public Chebyshev(Vector template, double eigmin, double eigmax) {
        p = template.copy();
        z = template.copy();
        r = template.copy();
        q = template.copy();
        setEigenvalues(eigmin, eigmax);
    }

    /**
     * Sets the eigenvalue estimates.
     * 
     * @param eigmin
     *            Smallest eigenvalue. Must be positive
     * @param eigmax
     *            Largest eigenvalue. Must be positive
     */
    public void setEigenvalues(double eigmin, double eigmax) {
        this.eigmin = eigmin;
        this.eigmax = eigmax;

        if (eigmin <= 0)
            throw new IllegalArgumentException("eigmin <= 0");
        if (eigmax <= 0)
            throw new IllegalArgumentException("eigmax <= 0");
        if (eigmin > eigmax)
            throw new IllegalArgumentException("eigmin > eigmax");
    }

	public Vector solve(Matrix A, Vector b, Vector x) throws IterativeSolverNotConvergedException {
       checkSizes(A, b, x);

       double alpha = 0, beta = 0, c = 0, d = 0;

       A.multAdd(-1, x, r.set(b));

       c = (eigmax - eigmin) / 2.0;
       d = (eigmax + eigmin) / 2.0;

       for (iter.setFirst(); !iter.converged(r, x); iter.next()) {
           M.apply(r, z);

           if (iter.isFirst()) {
               p.set(z);
               alpha = 2.0 / d;
           } else {
               beta = (alpha * c) / 2.0;
               beta *= beta;
               alpha = 1.0 / (d - beta);
               p.scale(beta).add(z);
           }

           A.mult(p, q);
           x.add(alpha, p);
           r.add(-alpha, q);
       }

       return x;
   }

}