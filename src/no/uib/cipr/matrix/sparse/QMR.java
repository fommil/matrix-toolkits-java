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
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;

/**
 * Quasi-Minimal Residual method. QMR solves the unsymmetric linear system
 * <code>Ax = b</code> using the Quasi-Minimal Residual method. QMR uses two
 * preconditioners, and by default these are the same preconditioner.
 * 
 * @author Templates
 */
public class QMR extends AbstractIterativeSolver {

    /**
     * Left preconditioner
     */
    private Preconditioner M1;

    /**
     * Right preconditioner
     */
    private Preconditioner M2;

    /**
     * Vectors for use in the iterative solution process
     */
    private Vector r, y, z, v, w, p, q, d, s, v_tld, w_tld, y_tld, z_tld,
            p_tld;

    /**
     * Constructor for QMR. Uses the given vector as template for creating
     * scratch vectors. Typically, the solution or the right hand side vector
     * can be passed, and the template is not modified
     * 
     * @param template
     *            Vector to use as template for the work vectors needed in the
     *            solution process
     */
    public QMR(Vector template) {
        M1 = M;
        M2 = M;
        r = template.copy();
        y = template.copy();
        z = template.copy();
        v = template.copy();
        w = template.copy();
        p = template.copy();
        q = template.copy();
        d = template.copy();
        s = template.copy();
        v_tld = template.copy();
        w_tld = template.copy();
        y_tld = template.copy();
        z_tld = template.copy();
        p_tld = template.copy();
    }

    /**
     * Constructor for QMR. Uses the given vector as template for creating
     * scratch vectors. Typically, the solution or the right hand side vector
     * can be passed, and the template is not modified. Allows setting different
     * right and left preconditioners
     * 
     * @param template
     *            Vector to use as template for the work vectors needed in the
     *            solution process
     * @param M1
     *            Left preconditioner
     * @param M2
     *            Right preconditioner
     */
    public QMR(Vector template, Preconditioner M1, Preconditioner M2) {
        this.M1 = M1;
        this.M2 = M2;
        r = template.copy();
        y = template.copy();
        z = template.copy();
        v = template.copy();
        w = template.copy();
        p = template.copy();
        q = template.copy();
        d = template.copy();
        s = template.copy();
        v_tld = template.copy();
        w_tld = template.copy();
        y_tld = template.copy();
        z_tld = template.copy();
        p_tld = template.copy();
    }

    public Vector solve(Matrix A, Vector b, Vector x)
            throws IterativeSolverNotConvergedException {
        checkSizes(A, b, x);

        double rho = 0, rho_1 = 0, xi = 0, gamma = 1., gamma_1 = 0, theta = 0, theta_1 = 0, eta = -1., delta = 0, ep = 0, beta = 0;

        A.multAdd(-1, x, r.set(b));

        v_tld.set(r);
        M1.apply(v_tld, y);
        rho = y.norm(Norm.Two);

        w_tld.set(r);
        M2.transApply(w_tld, z);
        xi = z.norm(Norm.Two);

        for (iter.setFirst(); !iter.converged(r, x); iter.next()) {

            if (rho == 0)
                throw new IterativeSolverNotConvergedException(
                        NotConvergedException.Reason.Breakdown, "rho", iter);

            if (xi == 0)
                throw new IterativeSolverNotConvergedException(
                        NotConvergedException.Reason.Breakdown, "xi", iter);

            v.set(1 / rho, v_tld);
            y.scale(1 / rho);
            w.set(1 / xi, w_tld);
            z.scale(1 / xi);

            delta = z.dot(y);

            if (delta == 0)
                throw new IterativeSolverNotConvergedException(
                        NotConvergedException.Reason.Breakdown, "delta", iter);

            M2.apply(y, y_tld);
            M1.transApply(z, z_tld);

            if (iter.isFirst()) {
                p.set(y_tld);
                q.set(z_tld);
            } else {
                p.scale(-xi * delta / ep).add(y_tld);
                q.scale(-rho * delta / ep).add(z_tld);
            }

            A.mult(p, p_tld);

            ep = q.dot(p_tld);

            if (ep == 0)
                throw new IterativeSolverNotConvergedException(
                        NotConvergedException.Reason.Breakdown, "ep", iter);

            beta = ep / delta;

            if (beta == 0)
                throw new IterativeSolverNotConvergedException(
                        NotConvergedException.Reason.Breakdown, "beta", iter);

            v_tld.set(-beta, v).add(p_tld);
            M1.apply(v_tld, y);
            rho_1 = rho;
            rho = y.norm(Norm.Two);

            A.transMultAdd(q, w_tld.set(-beta, w));
            M2.transApply(w_tld, z);
            xi = z.norm(Norm.Two);

            gamma_1 = gamma;
            theta_1 = theta;
            theta = rho / (gamma_1 * beta);
            gamma = 1 / Math.sqrt(1 + theta * theta);

            if (gamma == 0)
                throw new IterativeSolverNotConvergedException(
                        NotConvergedException.Reason.Breakdown, "gamma", iter);

            eta = -eta * rho_1 * gamma * gamma / (beta * gamma_1 * gamma_1);

            if (iter.isFirst()) {
                d.set(eta, p);
                s.set(eta, p_tld);
            } else {
                double val = theta_1 * theta_1 * gamma * gamma;
                d.scale(val).add(eta, p);
                s.scale(val).add(eta, p_tld);
            }

            x.add(d);
            r.add(-1, s);
        }

        return x;
    }

    @Override
    public void setPreconditioner(Preconditioner M) {
        super.setPreconditioner(M);
        M1 = M;
        M2 = M;
    }

}
