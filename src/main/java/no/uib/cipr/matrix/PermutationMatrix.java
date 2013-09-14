package no.uib.cipr.matrix;

import com.github.fommil.netlib.LAPACK;

import java.util.BitSet;

/**
 * Matrix that represents a permutation of another matrix's
 * rows / columns.
 * <p>
 * NOTE: the transpose of a permutation matrix is its inverse.
 *
 * @author Sam Halliday
 */
public class PermutationMatrix extends AbstractMatrix {

  /**
   * The sequential row permutations to perform, e.g. (2, 3, 3) means:
   * permute row 1 with row 2, then permute row 2 with row 3,
   * then permute row 3 with row 3 (i.e. do nothing).
   * <p>
   * Using this factory will ensure that LAPACK optimisations are
   * available for multiplication operations.
   *
   * @param pivots using fortran (1-indexed) notation.
   */
  public static PermutationMatrix fromPartialPivots(int pivots[]) {
    int[] permutations = new int[pivots.length];
    for (int i = 0; i < pivots.length; i++) {
      permutations[i] = i;
    }

    for (int i = 0; i < pivots.length; i++) {
      int j = pivots[i] - 1;
      if (j == i) continue;
      int tmp = permutations[i];
      permutations[i] = permutations[j];
      permutations[j] = tmp;
    }

    return new PermutationMatrix(permutations, pivots);
  }

  private int[] permutations, pivots;

  private boolean transposed;

  // the instantaneous permutations to perform (zero-indexed)
  // http://en.wikipedia.org/wiki/Permutation_matrix
  public PermutationMatrix(int permutations[]) {
    this(permutations, null);
  }

  // permutations - instantaneous (zero-indexed)
  // pivots - sequential (fortran-indexed)
  private PermutationMatrix(int permutations[], int pivots[]) {
    super(permutations.length, permutations.length);
    this.permutations = permutations;
    BitSet bitset = new BitSet();
    for (int i : permutations) {
      if (bitset.get(i))
        throw new IllegalArgumentException("non-unique permutations: " + i);
      bitset.set(i);
    }
    this.pivots = pivots;
  }

  @Override
  public double get(int row, int column) {
    if (!transposed && permutations[row] == column) return 1;
    if (transposed && permutations[column] == row) return 1;
    return 0;
  }

  @Override
  public Matrix transpose() {
    transposed = !transposed;
    return this;
  }

  @Override
  public Matrix mult(Matrix B, Matrix C) {
    if (C instanceof DenseMatrix) return mult(B, (DenseMatrix) C);
    return super.mult(B, C);
  }

  public Matrix mult(Matrix B, DenseMatrix C) {
    if (pivots == null) return super.mult(B, C);
    checkMultAdd(B, C);
    C.set(B);

    LAPACK.getInstance().dlaswp(C.numColumns(), C.getData(), Matrices.ld(C.numRows()), 1, pivots.length, pivots, transposed? -1 : 1);
    return C;
  }

  @Override
  public Matrix transAmult(Matrix B, Matrix C) {
    if (C instanceof DenseMatrix) return transAmult(B, (DenseMatrix) C);
    return super.transAmult(B, C);
  }

  public Matrix transAmult(Matrix B, DenseMatrix C) {
    if (pivots == null) return super.transAmult(B, C);
    checkTransAmultAdd(B, C);
    C.set(B);

    LAPACK.getInstance().dlaswp(C.numColumns(), C.getData(), Matrices.ld(C.numRows()), 1, pivots.length, pivots, transposed? 1 : -1);
    return C;
  }

}
