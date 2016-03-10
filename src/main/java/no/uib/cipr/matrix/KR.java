/*
 * Copyright (C) 2015 Rogério Pontes
 */

package no.uib.cipr.matrix;

/**
 * Class to compute The Khatri Rao multiplicaiton given a Matrix A and B Khatri
 * Rao Multiplication <code>C = A*B</code>
 * 
 */
public class KR {

    /**
     * Left Matrix of the computation
     */
    private final Matrix A;

    /**
     * Right Matrix of the compution
     */
    private final Matrix B;

    /**
     * @param A
     * @param B
     *            Matrix such that <code>B.numColumns() == A.numColumns()</code>
     */
    public KR(Matrix A, Matrix B) {
        checkKhatriRaoArguments(A, B);
        this.A = A;
        this.B = B;
    }

    /**
     * Tests if the matrices have an equal number of columns.
     * 
     * @param A
     * @param B
     */
    private static void checkKhatriRaoArguments(Matrix A, Matrix B) {
        if (A.numColumns() != B.numColumns())
            throw new IndexOutOfBoundsException(
                    "A.numColumns != B.numColumns (" + A.numColumns() + " != "
                            + B.numColumns() + ")");
    }

    /**
     * @param C
     *            Matrix such that <code>C.numRows() == A.numRows()*B.numRows()
     *    </code> and <code>B.numColumns() == C.numColumns()</code> Checks the
     *            dimensions of the result of <code>khatriRao</code>
     */
    private void checkKhatriRao(Matrix C) {
        if (A.numColumns() != C.numColumns())
            throw new IndexOutOfBoundsException(
                    "A.numColumns != C.numColumns (" + A.numColumns() + " != "
                            + C.numColumns() + ")");
        if (A.numRows() * B.numRows() != C.numRows()) {
            throw new IndexOutOfBoundsException(
                    "C.numRows != A.numRows * B.numRows ( " + C.numRows()
                            + " != " + A.numRows() * B.numRows() + " )");
        }
    }

    /**
     * 
     * @param C
     *            The result of the Khatri Rao Multiplication of A*B
     * @return C
     */
    public Matrix multiply(Matrix C) {
        checkKhatriRao(C);

        for (int i = 0; i < A.numColumns(); ++i)
            for (int j = 0; j < A.numRows(); ++j)
                for (int k = 0; k < B.numRows(); ++k) {
                    double value = A.get(j, i) * B.get(k, i);
                    int destLine = B.numRows() * j + k;
                    C.add(destLine, i, value);
                }
        return C;

    }

}
