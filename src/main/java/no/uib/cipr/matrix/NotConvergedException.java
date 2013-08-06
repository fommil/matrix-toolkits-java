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

package no.uib.cipr.matrix;

/**
 * Signals lack of convergence of an iterative process
 */
public class NotConvergedException extends Exception {

    private static final long serialVersionUID = -2305369220010776320L;

    /**
     * Possible reasons for lack of convergence
     */
    public enum Reason {

        /**
         * Did not converge after a maximum number of iterations
         */
        Iterations,

        /**
         * Divergence detected
         */
        Divergence,

        /**
         * The iterative process detected a breakdown
         */
        Breakdown
    }

    /**
     * The reason for this exception
     */
    protected Reason reason;

    /**
     * Constructor for NotConvergedException
     * 
     * @param reason
     *            The reason for the lack of convergence
     * @param message
     *            A more descriptive explanation
     */
    public NotConvergedException(Reason reason, String message) {
        super(message);
        this.reason = reason;
    }

    /**
     * Constructor for NotConvergedException. No message is provided
     * 
     * @param reason
     *            The reason for the lack of convergence
     */
    public NotConvergedException(Reason reason) {
        this.reason = reason;
    }

    /**
     * Returns the reason for the exception
     */
    public Reason getReason() {
        return reason;
    }
}
