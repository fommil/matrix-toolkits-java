/*
 * Copyright (C) 2006 Rafael de Pelegrini Soares
 *
 * MTJ additions.
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

import junit.framework.TestCase;

/**
 * QRP test
 */
public class QRPTest extends TestCase {

	public void testIdentity() {
		Matrix A = Matrices.identity(5);

		QRP qrp = QRP.factorize(A);

		assertEquals(A, mult(qrp.getQ(), qrp.getPVector(), qrp.getR(), A.copy().zero()));
	}

	public void testRectangularIdentity1() {
		Matrix A = new DenseMatrix(5, 6);
		for (MatrixEntry i : A) {
			if (i.row() == i.column())
				i.set(1);
		}

		QRP qrp = QRP.factorize(A);

		assertEquals(A, mult(qrp.getQ(), qrp.getPVector(), qrp.getR(), A.copy().zero()));
	}

	public void testRectangularIdentity2() {
		Matrix A = new DenseMatrix(6, 5);
		for (MatrixEntry i : A) {
			if (i.row() == i.column())
				i.set(1);
		}

		QRP qrp = QRP.factorize(A);

		assertEquals(A, mult(qrp.getQ(), qrp.getPVector(), qrp.getR(), A.copy().zero()));
	}

	public void testOrthogonality() {
		Matrix A = Matrices.random(6, 4);

		QRP qrp = QRP.factorize(A);

		Matrix QP = qrp.getQ();

		assertEquals(Matrices.identity(QP.numRows()), QP.transAmult(QP, QP.copy().zero()));
	}

	public void testOrthogonality2() {
		Matrix A = Matrices.random(4, 6);

		QRP qrp = QRP.factorize(A);

		Matrix QP = qrp.getQ();

		assertEquals(Matrices.identity(QP.numRows()), QP.transAmult(QP, QP.copy().zero()));
	}

	public void testRandSquare() {
		Matrix A = Matrices.random(5, 5);

		QRP qrp = QRP.factorize(A);
		Matrix Q = qrp.getQ();
		Matrix R = qrp.getR();
		int P[] = qrp.getPVector();

		assertEquals(A, mult(Q, P, R, A.copy().zero()));
	}

	public void testRandRectangular1() {
		Matrix A = Matrices.random(4, 6);

		QRP qrp = QRP.factorize(A);

		assertEquals(A, mult(qrp.getQ(), qrp.getPVector(), qrp.getR(), A.copy().zero()));
	}

	public void testRandRectangular2() {
		Matrix A = Matrices.random(6, 4);

		QRP qrp = QRP.factorize(A);

		assertEquals(A, mult(qrp.getQ(), qrp.getPVector(), qrp.getR(), A.copy().zero()));
	}

	public void testPivotingMatrix() {
		Matrix A = Matrices.random(6, 4);

		QRP qrp = QRP.factorize(A);
		Matrix Q = qrp.getQ();
		Matrix R = qrp.getR();
		Matrix P = qrp.getP();

		Matrix C = A.copy().zero();
		Matrix D = A.copy().zero();

		Q.mult(R, C);
		C.transBmult(P, D);

		assertEquals(A, D);
	}

	public void testRank1() {
		Matrix rand = Matrices.random(6, 4);

		Matrix A = new DenseMatrix(rand.numRows(), rand.numRows());
		rand.transBmult(rand, A);

		QRP qrp = QRP.factorize(A);

		assertEquals(Math.min(rand.numRows(), rand.numColumns()), qrp.getRank());
	}

	public void testRank2() {
		Matrix rand = Matrices.random(4, 6);

		Matrix A = new DenseMatrix(rand.numRows(), rand.numRows());
		rand.transBmult(rand, A);

		QRP qrp = QRP.factorize(A);

		assertEquals(Math.min(rand.numRows(), rand.numColumns()), qrp.getRank());
	}

	/**
	 * Executes the multiplication C = Q*P*R
	 * @param Q matrix Q
	 * @param P column permutation of R
	 * @param R matrix R
	 * @param C matrix tu put the results
	 * @return the matrix C
	 */
	protected Matrix mult(Matrix Q, int[] P, Matrix R, Matrix C) {
		for (int i = 0; i < Q.numRows(); ++i) {
			for (int j = 0; j < C.numColumns(); ++j) {
				double dot = 0;
				for (int k = 0; k < Q.numColumns(); ++k) {
					dot += Q.get(i, k) * R.get(k, j);
				}
				C.add(i, P[j], dot);
			}
		}
		return C;
	}
	private static final double DELTA = 1e-12;

	/**
	 * Assert that the given matrices are identical
	 */
	protected void assertEquals(Matrix A, Matrix B) {
		assertEquals(A.numRows(), B.numRows());
		assertEquals(A.numColumns(), B.numColumns());
		for (int i = 0; i < A.numRows(); ++i) {
			for (int j = 0; j < A.numColumns(); ++j) {
				assertEquals(A.get(i, j), B.get(i, j), DELTA);
			}
		}
	}
}
