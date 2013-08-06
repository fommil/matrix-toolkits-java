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
 * Basic matrix interface. It holds <code>double</code>s in a rectangular 2D
 * array, and it is used alongside <code>Vector</code> in numerical
 * computations. Implementing classes decides on the actual storage.
 * 
 * <h4>Basic operations</h4>
 * <p>
 * Use <code>numRows</code> and <code>numColumns</code> to get the basic
 * size of a matrix. <code>get(int,int)</code> gets an element, and there are
 * corresponding <code>set(int,int,double)</code> and
 * <code>add(int,int,double)</code> methods as well. Note that matrix indices
 * are zero-based (typical for Java and C). This means that the row-indices
 * range from 0 to <code>numRows-1</code>, likewise for the columns. It is
 * legal to have <code>numRows</code> or <code>numColumns</code> equal zero.
 * </p>
 * <p>
 * Other basic operations are <code>zero</code> which zeros all the entries of
 * the matrix, which can be cheaper than either zeroing the matrix manually, or
 * creating a new matrix, and the operation <code>copy</code> which creates a
 * deep copy of the matrix. This copy has separate storage, but starts with the
 * same contents as the current matrix.
 * </p>
 * 
 * <h4>Iterators</h4>
 * <p>
 * The matrix interface extends <code>Iterable</code>, and the iterator
 * returns a <code>MatrixEntry</code> which contains current index and entry
 * value. Note that the iterator may skip non-zero entries. Using an iterator,
 * many simple and efficient algorithms can be created. The iterator also
 * permits changing values in the matrix, however only non-zero entries can be
 * changed.
 * </p>
 * 
 * <h4>Basic linear algebra</h4>
 * <p>
 * A large selection of basic linear algebra operations are available. To ensure
 * high efficiency, little or no internal memory allocation is done, and the
 * user is required to supply the output arguments.
 * </p>
 * <p>
 * The operations available include:
 * </p>
 * <dl>
 * <dt><i>Additions </i></dt>
 * <dd>Matrices can be added to each other, even if their underlying matrix
 * structures are different.</dd>
 * <dt><i>Multiplications </i></dt>
 * <dd>A matrix can be multiplied with vectors and other matrices. For
 * increased efficiency, a multiplication can be combined with addition and
 * scaling, and transpose matrix multiplications are also available.</dd>
 * <dt><i>Rank-updates </i></dt>
 * <dd>A matrix can be efficiently updated using low-rank updates. The updates
 * can be contained in both matrices or vectors.</dd>
 * <dt><i>Transpositions </i></dt>
 * <dd>In-place transpositions of square matrices is supported, and the
 * transpose of a matrix can be stored in another matrix of compatible size
 * (possibly non-rectangular)</dd>
 * <dt><i>Solvers </i></dt>
 * <dd>Many dense and structured sparse matrices have fast, direct solvers, and
 * can be used to solve linear systems without creating a factorization. These
 * solvers are typically backed by subroutines in LAPACK</dd>
 * </dl>
 */
public interface Matrix extends Iterable<MatrixEntry> {

    /**
     * Number of rows in the matrix
     */
    int numRows();

    /**
     * Number of columns in the matrix
     */
    int numColumns();

    /**
     * Returns true if the matrix is square
     */
    boolean isSquare();

    /**
     * <code>A(row,column) = value</code>
     */
    void set(int row, int column, double value);

    /**
     * <code>A(row,column) += value</code>
     */
    void add(int row, int column, double value);

    /**
     * Returns <code>A(row,column)</code>
     */
    double get(int row, int column);

    /**
     * Creates a deep copy of the matrix
     * 
     * @return A
     */
    Matrix copy();

    /**
     * Zeros all the entries in the matrix, while preserving any underlying
     * structure. Useful for general, unstructured matrices.
     * 
     * @return A
     */
    Matrix zero();

    /**
     * <code>y = A*x</code>
     * 
     * @param x
     *            Vector of size <code>A.numColumns()</code>
     * @param y
     *            Vector of size <code>A.numRows()</code>
     * @return y
     */
    Vector mult(Vector x, Vector y);

    /**
     * <code>y = alpha*A*x</code>
     * 
     * @param x
     *            Vector of size <code>A.numColumns()</code>
     * @param y
     *            Vector of size <code>A.numRows()</code>
     * @return y
     */
    Vector mult(double alpha, Vector x, Vector y);

    /**
     * <code>y = A*x + y</code>
     * 
     * @param x
     *            Vector of size <code>A.numColumns()</code>
     * @param y
     *            Vector of size <code>A.numRows()</code>
     * @return y
     */
    Vector multAdd(Vector x, Vector y);

    /**
     * <code>y = alpha*A*x + y</code>
     * 
     * @param x
     *            Vector of size <code>A.numColumns()</code>
     * @param y
     *            Vector of size <code>A.numRows()</code>
     * @return y
     */
    Vector multAdd(double alpha, Vector x, Vector y);

    /**
     * <code>y = A<sup>T</sup>*x</code>
     * 
     * @param x
     *            Vector of size <code>A.numRows()</code>
     * @param y
     *            Vector of size <code>A.numColumns()</code>
     * @return y
     */
    Vector transMult(Vector x, Vector y);

    /**
     * <code>y = alpha*A<sup>T</sup>*x</code>
     * 
     * @param x
     *            Vector of size <code>A.numRows()</code>
     * @param y
     *            Vector of size <code>A.numColumns()</code>
     * @return y
     */
    Vector transMult(double alpha, Vector x, Vector y);

    /**
     * <code>y = A<sup>T</sup>*x + y</code>
     * 
     * @param x
     *            Vector of size <code>A.numRows()</code>
     * @param y
     *            Vector of size <code>A.numColumns()</code>
     * @return y
     */
    Vector transMultAdd(Vector x, Vector y);

    /**
     * <code>y = alpha*A<sup>T</sup>*x + y</code>
     * 
     * @param x
     *            Vector of size <code>A.numRows()</code>
     * @param y
     *            Vector of size <code>A.numColumns()</code>
     * @return y
     */
    Vector transMultAdd(double alpha, Vector x, Vector y);

    /**
     * <code>x = A\b</code>. Not all matrices support this operation, those
     * that do not throw <code>UnsupportedOperationException</code>. Note
     * that it is often more efficient to use a matrix decomposition and its
     * associated solver
     * 
     * @param b
     *            Vector of size <code>A.numRows()</code>
     * @param x
     *            Vector of size <code>A.numColumns()</code>
     * @return x
     * @throws MatrixSingularException
     *             If the matrix is singular
     * @throws MatrixNotSPDException
     *             If the solver assumes that the matrix is symmetrical,
     *             positive definite, but that that property does not hold
     */
    Vector solve(Vector b, Vector x) throws MatrixSingularException,
            MatrixNotSPDException;

    /**
     * <code>x = A<sup>T</sup>\b</code>. Not all matrices support this
     * operation, those that do not throw
     * <code>UnsupportedOperationException</code>. Note that it is often more
     * efficient to use a matrix decomposition and its associated solver
     * 
     * @param b
     *            Vector of size <code>A.numColumns()</code>
     * @param x
     *            Vector of size <code>A.numRows()</code>
     * @return x
     * @throws MatrixSingularException
     *             If the matrix is singular
     * @throws MatrixNotSPDException
     *             If the solver assumes that the matrix is symmetrical,
     *             positive definite, but that that property does not hold
     */
    Vector transSolve(Vector b, Vector x) throws MatrixSingularException,
            MatrixNotSPDException;

    /**
     * <code>A = x*x<sup>T</sup> + A</code>. The matrix must be square, and
     * the vector of the same length
     * 
     * @return A
     */
    Matrix rank1(Vector x);

    /**
     * <code>A = alpha*x*x<sup>T</sup> + A</code>. The matrix must be
     * square, and the vector of the same length
     * 
     * @return A
     */
    Matrix rank1(double alpha, Vector x);

    /**
     * <code>A = x*y<sup>T</sup> + A</code>. The matrix must be square, and
     * the vectors of the same length
     * 
     * @return A
     */
    Matrix rank1(Vector x, Vector y);

    /**
     * <code>A = alpha*x*y<sup>T</sup> + A</code>. The matrix must be
     * square, and the vectors of the same length
     * 
     * @return A
     */
    Matrix rank1(double alpha, Vector x, Vector y);

    /**
     * <code>A = x*y<sup>T</sup> + y*x<sup>T</sup> + A</code>. The matrix
     * must be square, and the vectors of the same length
     * 
     * @return A
     */
    Matrix rank2(Vector x, Vector y);

    /**
     * <code>A = alpha*x*y<sup>T</sup> + alpha*y*x<sup>T</sup> + A</code>.
     * The matrix must be square, and the vectors of the same length
     * 
     * @return A
     */
    Matrix rank2(double alpha, Vector x, Vector y);

    /**
     * <code>C = A*B</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix mult(Matrix B, Matrix C);

    /**
     * <code>C = alpha*A*B</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix mult(double alpha, Matrix B, Matrix C);

    /**
     * <code>C = A*B + C</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix multAdd(Matrix B, Matrix C);

    /**
     * <code>C = alpha*A*B + C</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix multAdd(double alpha, Matrix B, Matrix C);

    /**
     * <code>C = A<sup>T</sup>*B</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix transAmult(Matrix B, Matrix C);

    /**
     * <code>C = alpha*A<sup>T</sup>*B</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix transAmult(double alpha, Matrix B, Matrix C);

    /**
     * <code>C = A<sup>T</sup>*B + C</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix transAmultAdd(Matrix B, Matrix C);

    /**
     * <code>C = alpha*A<sup>T</sup>*B + C</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix transAmultAdd(double alpha, Matrix B, Matrix C);

    /**
     * <code>C = A*B<sup>T</sup></code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix transBmult(Matrix B, Matrix C);

    /**
     * <code>C = alpha*A*B<sup>T</sup></code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix transBmult(double alpha, Matrix B, Matrix C);

    /**
     * <code>C = A*B<sup>T</sup> + C</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix transBmultAdd(Matrix B, Matrix C);

    /**
     * <code>C = alpha*A*B<sup>T</sup> + C</code>
     * 
     * @param B
     *            Matrix such that <code>B.numRows() == A.numRows()</code> and
     *            <code>B.numColumns() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numColumns() == C.numColumns()</code>
     * @return C
     */
    Matrix transBmultAdd(double alpha, Matrix B, Matrix C);

    /**
     * <code>C = A<sup>T</sup>*B<sup>T</sup></code>
     * 
     * @param B
     *            Matrix such that <code>B.numColumns() == A.numRows()</code>
     *            and <code>B.numRows() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numRows() == C.numColumns()</code>
     * @return C
     */
    Matrix transABmult(Matrix B, Matrix C);

    /**
     * <code>C = alpha*A<sup>T</sup>*B<sup>T</sup></code>
     * 
     * @param B
     *            Matrix such that <code>B.numColumns() == A.numRows()</code>
     *            and <code>B.numRows() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numRows() == C.numColumns()</code>
     * @return C
     */
    Matrix transABmult(double alpha, Matrix B, Matrix C);

    /**
     * <code>C = A<sup>T</sup>*B<sup>T</sup> + C</code>
     * 
     * @param B
     *            Matrix such that <code>B.numColumns() == A.numRows()</code>
     *            and <code>B.numRows() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numRows() == C.numColumns()</code>
     * @return C
     */
    Matrix transABmultAdd(Matrix B, Matrix C);

    /**
     * <code>C = alpha*A<sup>T</sup>*B<sup>T</sup> + C</code>
     * 
     * @param B
     *            Matrix such that <code>B.numColumns() == A.numRows()</code>
     *            and <code>B.numRows() == C.numColumns()</code>
     * @param C
     *            Matrix such that <code>C.numRows() == A.numColumns()</code>
     *            and <code>B.numRows() == C.numColumns()</code>
     * @return C
     */
    Matrix transABmultAdd(double alpha, Matrix B, Matrix C);

    /**
     * <code>X = A\B</code>. Not all matrices support this operation, those
     * that do not throw <code>UnsupportedOperationException</code>. Note
     * that it is often more efficient to use a matrix decomposition and its
     * associated solver
     * 
     * @param B
     *            Matrix with the same number of rows as <code>A</code>, and
     *            the same number of columns as <code>X</code>
     * @param X
     *            Matrix with a number of rows equal <code>A.numColumns()</code>,
     *            and the same number of columns as <code>B</code>
     * @return X
     * @throws MatrixSingularException
     *             If the matrix is singular
     * @throws MatrixNotSPDException
     *             If the solver assumes that the matrix is symmetrical,
     *             positive definite, but that that property does not hold
     */
    Matrix solve(Matrix B, Matrix X) throws MatrixSingularException,
            MatrixNotSPDException;

    /**
     * <code>X = A<sup>T</sup>\B</code>. Not all matrices support this
     * operation, those that do not throw
     * <code>UnsupportedOperationException</code>. Note that it is often more
     * efficient to use a matrix decomposition and its associated transpose
     * solver
     * 
     * @param B
     *            Matrix with a number of rows equal <code>A.numColumns()</code>,
     *            and the same number of columns as <code>X</code>
     * @param X
     *            Matrix with the same number of rows as <code>A</code>, and
     *            the same number of columns as <code>B</code>
     * @return X
     * @throws MatrixSingularException
     *             If the matrix is singular
     * @throws MatrixNotSPDException
     *             If the solver assumes that the matrix is symmetrical,
     *             positive definite, but that that property does not hold
     */
    Matrix transSolve(Matrix B, Matrix X) throws MatrixSingularException,
            MatrixNotSPDException;

    /**
     * <code>A = C*C<sup>T</sup> + A</code>. The matrices must be square
     * and of the same size
     * 
     * @return A
     */
    Matrix rank1(Matrix C);

    /**
     * <code>A = alpha*C*C<sup>T</sup> + A</code>. The matrices must be
     * square and of the same size
     * 
     * @return A
     */
    Matrix rank1(double alpha, Matrix C);

    /**
     * <code>A = C<sup>T</sup>*C + A</code> The matrices must be square and
     * of the same size
     * 
     * @return A
     */
    Matrix transRank1(Matrix C);

    /**
     * <code>A = alpha*C<sup>T</sup>*C + A</code> The matrices must be
     * square and of the same size
     * 
     * @return A
     */
    Matrix transRank1(double alpha, Matrix C);

    /**
     * <code>A = B*C<sup>T</sup> + C*B<sup>T</sup> + A</code>. This
     * matrix must be square
     * 
     * @param B
     *            Matrix with the same number of rows as <code>A</code> and
     *            the same number of columns as <code>C</code>
     * @param C
     *            Matrix with the same number of rows as <code>A</code> and
     *            the same number of columns as <code>B</code>
     * @return A
     */
    Matrix rank2(Matrix B, Matrix C);

    /**
     * <code>A = alpha*B*C<sup>T</sup> + alpha*C*B<sup>T</sup> + A</code>.
     * This matrix must be square
     * 
     * @param B
     *            Matrix with the same number of rows as <code>A</code> and
     *            the same number of columns as <code>C</code>
     * @param C
     *            Matrix with the same number of rows as <code>A</code> and
     *            the same number of columns as <code>B</code>
     * @return A
     */
    Matrix rank2(double alpha, Matrix B, Matrix C);

    /**
     * <code>A = B<sup>T</sup>*C + C<sup>T</sup>*B + A</code>. This
     * matrix must be square
     * 
     * @param B
     *            Matrix with the same number of rows as <code>C</code> and
     *            the same number of columns as <code>A</code>
     * @param C
     *            Matrix with the same number of rows as <code>B</code> and
     *            the same number of columns as <code>A</code>
     * @return A
     */
    Matrix transRank2(Matrix B, Matrix C);

    /**
     * <code>A = alpha*B<sup>T</sup>*C + alpha*C<sup>T</sup>*B + A</code>.
     * This matrix must be square
     * 
     * @param B
     *            Matrix with the same number of rows as <code>C</code> and
     *            the same number of columns as <code>A</code>
     * @param C
     *            Matrix with the same number of rows as <code>B</code> and
     *            the same number of columns as <code>A</code>
     * @return A
     */
    Matrix transRank2(double alpha, Matrix B, Matrix C);

    /**
     * <code>A = alpha*A</code>
     * 
     * @return A
     */
    Matrix scale(double alpha);

    /**
     * <code>A=B</code>. The matrices must be of the same size
     * 
     * @return A
     */
    Matrix set(Matrix B);

    /**
     * <code>A=alpha*B</code>. The matrices must be of the same size
     * 
     * @return A
     */
    Matrix set(double alpha, Matrix B);

    /**
     * <code>A = B + A</code>. The matrices must be of the same size
     * 
     * @return A
     */
    Matrix add(Matrix B);

    /**
     * <code>A = alpha*B + A</code>. The matrices must be of the same size
     * 
     * @return A
     */
    Matrix add(double alpha, Matrix B);

    /**
     * Transposes the matrix in-place. In most cases, the matrix must be square
     * for this to work.
     * 
     * @return This matrix
     */
    Matrix transpose();

    /**
     * Sets the tranpose of this matrix into <code>B</code>. Matrix
     * dimensions must be compatible
     * 
     * @param B
     *            Matrix with as many rows as this matrix has columns, and as
     *            many columns as this matrix has rows
     * @return The matrix <code>B=A<sup>T</sup></code>
     */
    Matrix transpose(Matrix B);

    /**
     * Computes the given norm of the matrix
     * 
     * @param type
     *            The type of norm to compute
     */
    double norm(Norm type);

    /**
     * Supported matrix-norms. Note that <code>Maxvalue</code> is not a proper
     * matrix norm
     */
    enum Norm {

        /**
         * Maximum absolute row sum
         */
        One,

        /**
         * The root of sum of the sum of squares
         */
        Frobenius,

        /**
         * Maximum column sum
         */
        Infinity,

        /**
         * Largest entry in absolute value. Not a proper matrix norm
         */
        Maxvalue;

		/**
		 * @return the String as required by the netlib libraries to represent this norm.
		 */
		public String netlib() {
			// TODO: this is a bit of a hack
			// shouldn't need to know about the internals of netlib
		    if (this == One)
		        return "1";
		    else if (this == Infinity)
		        return "I";
		    else
		        throw new IllegalArgumentException(
		                "Norm must be the 1 or the Infinity norm");
		}

    }

}
