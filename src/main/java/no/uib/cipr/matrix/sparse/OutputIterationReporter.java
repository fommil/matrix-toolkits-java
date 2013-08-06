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

import java.io.OutputStream;
import java.io.PrintWriter;

import no.uib.cipr.matrix.Vector;

/**
 * Outputs iteration information to an output stream.
 */
public class OutputIterationReporter implements IterationReporter {

    /**
     * Platform-dependent output
     */
    private PrintWriter out;

    /**
     * Constructor for OutputIterationReporter
     * 
     * @param out
     *            Writes iteration count and current residual here
     */
    public OutputIterationReporter(OutputStream out) {
        this.out = new PrintWriter(out, true);
    }

    /**
     * Constructor for OutputIterationReporter, using <code>System.err</code>.
     */
    public OutputIterationReporter() {
        this(System.err);
    }

    public void monitor(double r, int i) {
        out.format("%10d % .12e\n", i, r);
        out.flush();
    }

    public void monitor(double r, Vector x, int i) {
        monitor(r, i);
    }

}
