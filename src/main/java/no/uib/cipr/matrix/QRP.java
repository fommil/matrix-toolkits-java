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

import com.github.fommil.netlib.LAPACK;
import org.netlib.util.intW;

/**
 * Computes QR decompositions with column pivoting:
 * 
 * {@code A*P = Q*R} where
 * 
 * {@code A(m,n)}, {@code Q(m,m)}, and {@code R(m,n)}, more generally:
 * 
 * {@code A*P = [Q1 Q2] * [R11, R12; 0 R22]} and {@code R22} elements are
 * negligible.
 * 
 */
public class QRP {

	/** Pivoting vector */
	int jpvt[];
	/**
	 * Scales for the reflectors
	 */
	final double[] tau;
	/**
	 * Factorisation sizes
	 */
	final int m,  n,  k;
	/** The factored matrix rank */
	int rank;
	/**
	 * Work array
	 */
	double[] work;
	/**
	 * The factored matrix
	 */
	final DenseMatrix Afact;
	/**
	 * The orthogonal matrix
	 */
	final DenseMatrix Q;
	/**
	 * The general upper triangular matrix.
	 */
	final DenseMatrix R;

	/**
	 * Constructs an empty QR decomposition
	 *
	 * @param m the number of rows.
	 * @param n the number of columns.
	 */
	public QRP(int m, int n) {
		this.m = m;
		this.n = n;
		this.k = Math.min(m, n);
		this.rank = 0;
		jpvt = new int[n];
		tau = new double[k];

		Q = new DenseMatrix(m, m);
		R = new DenseMatrix(m, n);
		Afact = new DenseMatrix(m, Math.max(m, n));

		int lwork1, lwork2;
		intW info = new intW(0);
		double dummy[] = new double[1];
		double ret[] = new double[1];

		LAPACK lapack = LAPACK.getInstance();

		// Query optimal workspace. First for computing the factorization
		lapack.dgeqrf(m, n, dummy, Matrices.ld(m), dummy, ret, -1, info);
		lwork1 = (info.val != 0) ? n : (int) ret[0];

		// Workspace needed for generating an explicit orthogonal matrix
		lapack.dorgqr(m, m, k, dummy, Matrices.ld(m), dummy, ret, -1, info);
		lwork2 = (info.val != 0) ? n : (int) ret[0];

		work = new double[Math.max(lwork1, lwork2)];
	}

	/**
	 * Convenience method to compute a QR decomposition
	 *
	 * @param A the matrix to decompose (not modified)
	 * @return Newly allocated decomposition
	 */
	public static QRP factorize(Matrix A) {
		return new QRP(A.numRows(), A.numColumns()).factor(A);
	}

	/**
	 * Executes a QR factorization for the given matrix.
	 *
	 * @param A the matrix to be factored (not modified)
	 * @return the factorization object
	 */
	public QRP factor(Matrix A) {
		if (Q.numRows() != A.numRows())
			throw new IllegalArgumentException("Q.numRows() != A.numRows()");
		else if (R.numColumns() != A.numColumns())
			throw new IllegalArgumentException("R.numColumns() != A.numColumns()");

		// copy A values in Afact
		Afact.zero();
		for (MatrixEntry e : A) {
			Afact.set(e.row(), e.column(), e.get());
		}

		intW info = new intW(0);
		LAPACK lapack = LAPACK.getInstance();

		/*
		 * Calculate factorisation
		 */
		lapack.dgeqp3(m, n, Afact.getData(), Matrices.ld(m), jpvt, tau, work,
			work.length, info);

		if (info.val < 0)
			throw new IllegalArgumentException();

		/*
		 * Get R from Afact
		 */
		R.zero();
		for (MatrixEntry e : Afact) {
			if (e.row() <= e.column() && e.column() < R.numColumns()) {
				R.set(e.row(), e.column(), e.get());
			}
		}

		/*
		 * Calculate the rank based on a precision EPS
		 */
		final double EPS = 1e-12;
		for (rank = 0; rank < k; rank++) {
			if (Math.abs(R.get(rank, rank)) < EPS)
				break;
		}

		/*
		 * Explicit the orthogonal matrix
		 */
		lapack.dorgqr(m, m, k, Afact.getData(), Matrices.ld(m), tau, work,
			work.length, info);
		for (MatrixEntry e : Afact) {
			if (e.column() < Q.numColumns())
				Q.set(e.row(), e.column(), e.get());
		}

		if (info.val < 0)
			throw new IllegalArgumentException();

		// Adjust the permutation to zero offset
		for (int i = 0; i < jpvt.length; i++) {
			--jpvt[i];
		}

		return this;
	}

	/**
	 * Returns the upper triangular factor
	 */
	public DenseMatrix getR() {
		return R;
	}

	/**
	 * Returns the orthogonal matrix
	 */
	public DenseMatrix getQ() {
		return Q;
	}

	/**
	 * Returns the column pivoting vector.
	 * This function is cheaper than {@link #getP()}.
	 */
	public int[] getPVector() {
		return jpvt;
	}

	/**
	 * Returns the column pivoting matrix.
	 * This function allocates a new Matrix to be returned,
	 * a more cheap option is tu use {@link #getPVector()}.
	 */
	public Matrix getP() {
		Matrix P = new DenseMatrix(jpvt.length, jpvt.length);
		for (int i = 0; i < jpvt.length; i++) {
			P.set(jpvt[i], i, 1);
		}
		return P;
	}

	/**
	 * Returns the rank of the factored matrix
	 */
	public int getRank() {
		return rank;
	}
}
