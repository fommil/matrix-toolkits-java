package no.uib.cipr.matrix;

import java.util.BitSet;

/**
 * Matrix that represents a permutation of another matrix's
 * rows / columns.
 *
 * @author Sam Halliday
 */
public class PermutationMatrix extends AbstractMatrix {

  /**
   * The sequential row permutations to perform, e.g. (2, 3, 3) means:
   * permute row 1 with row 2, then permute row 2 with row 3,
   * then permute row 3 with row 3 (i.e. do nothing).
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

  private final int[] permutations, pivots;

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
    if (permutations[row] == column) return 1;
    return 0;
  }

  // TODO: dlaswp for DenseMatrix and non-null pivots
  // LAPACK.getInstance().dlaswp(out.numColumns(), out.getData(), out.numRows(), 1, piv.length, piv, 1);

}
