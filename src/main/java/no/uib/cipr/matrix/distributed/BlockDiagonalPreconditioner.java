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

package no.uib.cipr.matrix.distributed;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.Preconditioner;

/**
 * Block diagonal preconditioner
 *
 * @deprecated the <code>no.uib.cipr.matrix.distributed</code> package has been deprecated because
 * of a number of hard to fix concurrency bugs. It is distributed only for backwards compatibility,
 * but is not recommended. The utility of this package is questionable, as it does not allow
 * distribution of computation between JVMs or across a network. For many people, distributed
 * computing of multiple matrices can be achieved at a user-level through the
 * <a href="http://jppf.org">JPPF Framework</a>.
 * Users who need to deal with few very large matrices may wish to implement their own storage classes
 * and solvers using JPPF, but this will not be supported directly in matrix-toolkits-java.
 */
@Deprecated
public class BlockDiagonalPreconditioner implements Preconditioner {

    /**
     * Preconditioner to apply on the block
     */
    private Preconditioner prec;

    /**
     * Constructor for BlockDiagonalPreconditioner
     * 
     * @param prec
     *            Preconditioner to apply on the blocks. As this preconditioner
     *            is meant to be used for distributed memory matrices, the
     *            preconditioner should be constructed on
     *            <code>A.getBlock()</code>
     */
    public BlockDiagonalPreconditioner(Preconditioner prec) {
        this.prec = prec;
    }

    public Vector apply(Vector b, Vector x) {
        if (!(b instanceof DistVector) || !(x instanceof DistVector))
            throw new IllegalArgumentException("Vectors must be DistVectors");

        return prec.apply(((DistVector) b).getLocal(), ((DistVector) x)
                .getLocal());
    }

    public Vector transApply(Vector b, Vector x) {
        if (!(b instanceof DistVector) || !(x instanceof DistVector))
            throw new IllegalArgumentException("Vectors must be DistVectors");

        return prec.transApply(((DistVector) b).getLocal(), ((DistVector) x)
                .getLocal());
    }

    public void setMatrix(Matrix A) {
        if (A instanceof DistRowMatrix)
            prec.setMatrix(((DistRowMatrix) A).getBlock());
        else if (A instanceof DistColMatrix)
            prec.setMatrix(((DistColMatrix) A).getBlock());
        else
            throw new IllegalArgumentException(
                    "!(A instanceof DistRowMatrix) && !(A instanceof DistColMatrix)");
    }
}
