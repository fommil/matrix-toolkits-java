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

import no.uib.cipr.matrix.DenseMatrix;

/**
 * Test of a dense matrix
 */
public class DenseMatrixTest extends MatrixTestAbstract {

    public DenseMatrixTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void createPrimary() throws Exception {
        int n = Utilities.getInt(1, max);
        int m = Utilities.getInt(1, max);
        A = new DenseMatrix(n, m);
        Ad = Utilities.populate(A);
    }

    @Override
    public void testMatrixSolve() {
        if (A.isSquare())
            super.testMatrixSolve();
    }

    @Override
    public void testTransMatrixSolve() {
        if (A.isSquare())
            super.testTransMatrixSolve();
    }

    @Override
    public void testTransVectorSolve() {
        if (A.isSquare())
            super.testTransVectorSolve();
    }

    @Override
    public void testVectorSolve() {
        if (A.isSquare())
            super.testVectorSolve();
    }

	public void testIssue13(){
		Vector bv = Matrices.random(100);
		Matrix am = Matrices.random(100, 50);
		Vector xv = new DenseVector(am.numColumns());
		for (int x = 0; x < am.numColumns(); x++) {
			xv.set(x, 1);
		}
		xv = Matrices.random(xv.size());
		xv = am.solve(bv, xv);
	}

	public void testIssue32(){

        // The issue here is that we should not allow matrices with more than
        // Integer.MAX_VALUE entries.
        boolean exceptionThrown = false;
        try
        {
            Matrix m = new DenseMatrix(Integer.MAX_VALUE, 2);
        }
        catch (IllegalArgumentException e)
        {
            exceptionThrown = true;
        }
        finally
        {
            assertTrue(exceptionThrown);
        }

        exceptionThrown = false;
        try
        {
            Matrix m = new DenseMatrix(Integer.MAX_VALUE, 3);
        }
        catch (IllegalArgumentException e)
        {
            exceptionThrown = true;
        }
        finally
        {
            assertTrue(exceptionThrown);
        }

        exceptionThrown = false;
        try
        {
            Matrix m = new DenseMatrix(Integer.MAX_VALUE - 10, 3);
        }
        catch (IllegalArgumentException e)
        {
            exceptionThrown = true;
        }
        finally
        {
            assertTrue(exceptionThrown);
        }
    }
}
