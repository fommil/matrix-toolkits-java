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

import no.uib.cipr.matrix.DenseLU;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import junit.framework.TestCase;

/**
 * Tests the dense LU decomposition
 */
public class DenseLUTest extends TestCase {

	/**
	 * Matrix to decompose
	 */
	private DenseMatrix A;

	private DenseMatrix I;

	private final int max = 100;

	public DenseLUTest(String arg0) {
		super(arg0);
	}

	@Override
	protected void setUp() throws Exception {
		int n = Utilities.getInt(1, max);
		A = new DenseMatrix(n, n);
		Utilities.populate(A);
		Utilities.addDiagonal(A, 1);
		while (Utilities.singular(A))
			Utilities.addDiagonal(A, 1);

		I = Matrices.identity(n);
	}

	@Override
	protected void tearDown() throws Exception {
		A = null;
		I = null;
	}

	public void testDenseLU() {
		int n = A.numRows();
		DenseLU lu = new DenseLU(n, n);
		lu.factor(A.copy());

		lu.solve(I);

		Matrix J = I.mult(A, new DenseMatrix(n, n));
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				if (i != j)
					assertEquals(J.get(i, j), 0, 1e-10);
				else
					assertEquals(J.get(i, j), 1, 1e-10);
	}

	public void testDenseLUtranspose() {
		int n = A.numRows();
		DenseLU lu = new DenseLU(n, n);
		lu.factor(A.copy());

		lu.transSolve(I);

		Matrix J = I.transAmult(A, new DenseMatrix(n, n));
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				if (i != j)
					assertEquals(J.get(i, j), 0, 1e-10);
				else
					assertEquals(J.get(i, j), 1, 1e-10);
	}

	public void testDenseLUrcond() {
		int n = A.numRows();
		DenseLU lu = new DenseLU(n, n);
		lu.factor(A.copy());

		lu.rcond(A, Matrix.Norm.One);
		lu.rcond(A, Matrix.Norm.Infinity);
	}

	public void testDenseLUToInput() {
		// MTJ bug in DenseLU code

		Matrix m = new DenseMatrix(3, 3);

		// -2 2 -3
		// -1 1 3
		// 2 0 -1
		m.set(0, 0, -2);
		m.set(0, 1, 2);
		m.set(0, 2, -3);
		m.set(1, 0, -1);
		m.set(1, 1, 1);
		m.set(1, 2, 3);
		m.set(2, 0, 2);
		m.set(2, 1, 0);
		m.set(2, 2, -1);

		// SHOULD BE:
		// L:
		// 1.000 0.000 0.000
		// -1.000 1.000 0.000
		// 0.500 0.000 1.000
		//
		// U:
		// -2.000 2.000 -3.000
		// 0.000 2.000 -4.000
		// 0.000 0.000 4.500
		//
		// Permutation matrix:
		// 1.000 0.000 0.000
		// 0.000 0.000 1.000
		// 0.000 1.000 0.000

		DenseLU dlu = DenseLU.factorize(m);

		// check that m = L . U
		Matrix lTimesU = new DenseMatrix(3, 3);
		dlu.getL().mult(dlu.getU(), lTimesU);
		int[] pivots = dlu.getPivots();
		for (MatrixEntry entry : m) {
			int row = entry.row();
			int col = entry.column();
			double val = entry.get();
			double valLU = pivots[row] * lTimesU.get(row, col);
			assert val == valLU : "Row " + row + ", Col " + col
					+ " wasn't equal! " + val + " " + valLU;
		}

		Matrix lu = dlu.getLU();
		// m == lu
		for (MatrixEntry entry : m) {
			int row = entry.row();
			int col = entry.column();
			double val = entry.get();
			double valLU = lu.get(row, col);
			assert val == valLU : "Row " + row + ", Col " + col
					+ " wasn't equal! " + val + " " + valLU;
		}
	}

}
