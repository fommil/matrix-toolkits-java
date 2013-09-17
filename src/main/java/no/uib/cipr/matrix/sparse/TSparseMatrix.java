package no.uib.cipr.matrix.sparse;

import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.io.MatrixVectorReader;

import java.io.IOException;

/**
 * GNU Trove backed general sparse matrix storage.
 * <p>
 * A good fit for unstructured sparse matrices, all
 * entries are stored using a hashmap constructed
 * from the row/column indices making mutations relatively
 * cheap
 * <p>
 * The trove hashes are constructed in such a way that
 * the iterator will go through rows (with potential
 * collisions on columns), giving
 * a slight performance optimisation when multiplying
 * other matrices or vectors (but at a significant cost
 * to transpose multiplications!)
 *
 * TODO: straight or transpose mult optimisation should be a user option... potentially both! (at memory cost)
 *
 * @author Sam Halliday
 */
public class TSparseMatrix extends AbstractMatrix {

  public TSparseMatrix(int numRows, int numColumns) {
    super(numRows, numColumns);
  }

  public TSparseMatrix(Matrix A) {
    super(A);
    // TODO
  }

  public TSparseMatrix(MatrixVectorReader r) throws IOException {
    super(0, 0);
    // TODO
  }
}