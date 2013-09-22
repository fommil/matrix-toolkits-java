package no.uib.cipr.matrix.sparse;

import lombok.AllArgsConstructor;
import lombok.ToString;
import lombok.extern.java.Log;
import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.io.MatrixInfo;
import no.uib.cipr.matrix.io.MatrixSize;
import no.uib.cipr.matrix.io.MatrixVectorReader;

import java.io.IOException;
import java.util.Iterator;

/**
 * A Linked List (with shortcuts to important nodes)
 * implementation of an {@code n x m} Matrix with {@code z} elements that
 * has a typical {@code O(z / m)} insertion / lookup cost
 * and an iterator that traverses columns then rows:
 * a good fit for unstructured sparse matrices.
 * <p/>
 * However, transpose
 * multiplications are costly and memory requirements
 * ({@code 2 int, 1 ref, 1 double} per matrix element, plus {@code n ref}s
 * for the cache) are higher than structured sparse matrix storage.
 * <p/>
 * TODO: double linkage option to remove transpose constraint
 *
 * @author Sam Halliday
 */
@Log
public class LinkedSparseMatrix extends AbstractMatrix {

  // java.util.LinkedList is doubly linked and therefore too heavyweight.
  @AllArgsConstructor
  @ToString(exclude = "tail")
  private static class Node {
    private final int row, col;
    private double val;
    private Node tail;
  }

  private class Linked {

    Node node;

    Node[] rows = new Node[numRows];

    public double get(int row, int col) {
      if (node != null && node.row == row && node.col == col) return node.val;
      Node prec = findPreceeding(row, col);
      if (prec == null) return 0;
      if (prec.tail == null) return 0;
      if (prec.tail.row == row && prec.tail.col == col) return prec.tail.val;
      return 0;
    }

    public void set(int row, int col, double val) {
      if (val == 0) {
        remove(row, col);
        return;
      }
      Node prec = findPreceeding(row, col);
      assert prec == null || prec.row < row || (prec.row == row && prec.col < col) : row + " " + col + " " + prec;
      if (node != null && node.row == row && node.col == col)
        node.val = val;
      else if (prec == null) {
        assert node == null || node.row > row || (node.row == row && node.col > col) : row + " " + col + " " + node;
        node = new Node(row, col, val, node);
      } else if (prec.tail == null || prec.tail.row != row || prec.tail.col != col) {
        prec.tail = new Node(row, col, val, prec.tail);
      } else {
        assert prec.tail.row == row && prec.tail.col == col;
        prec.tail.val = val;
      }
      invalidate(row, col);
    }

    public void remove(int row, int col) {
      Node prec = findPreceeding(row, col);
      if (prec == null) return;
      if (prec.tail == null) return;
      if (prec.tail.row != row || prec.tail.col != col) return;
      prec.tail = prec.tail.tail;
      invalidate(row, col);
    }

    void invalidate(int row, int col) {
      rows[row] = null;
    }

    // contains a pointer to the *last* Node for the row
    Node cached(int row) {
      for (int i = row; i > 0; i--) {
        Node cached = rows[row - 1];
        if (cached != null) {
          assert cached.row < row;
          return cached;
        }
      }
      return null;
    }

    // returns the node that either references this
    // index, or would be the one to update should it
    // be inserted: null if there are no entries.
    Node findPreceeding(int row, int col) {
      assert row >= 0 && col >= 0;
      Node last = row > 0 ? cached(row - 1) : null;
      Node cur = last != null ? last : node;
      while (cur != null && cur.row <= row) {
        if (cur.row == row && cur.col >= col) return last;
        last = cur;
        cur = cur.tail;
        if (cur != null && cur.row > last.row)
          rows[last.row] = last;
      }
      return last;
    }

  }

  private Linked rows;

  protected LinkedSparseMatrix(int numRows, int numColumns) {
    super(numRows, numColumns);
    rows = new Linked();
  }

  protected LinkedSparseMatrix(Matrix A) {
    super(A);
    rows = new Linked();
    set(A);
  }

  public LinkedSparseMatrix(MatrixVectorReader r) throws IOException {
    super(0, 0);
    try {
      MatrixInfo info = r.readMatrixInfo();
      if (info.isComplex()) throw new IllegalArgumentException("complex matrices not supported");
      if (!info.isCoordinate()) throw new IllegalArgumentException("only coordinate matrices supported");
      MatrixSize size = r.readMatrixSize(info);
      numRows = size.numRows();
      numColumns = size.numColumns();
      rows = new Linked();

      int nz = size.numEntries();
      int[] row = new int[nz];
      int[] column = new int[nz];
      double[] entry = new double[nz];
      r.readCoordinate(row, column, entry);
      r.add(-1, row);
      r.add(-1, column);
      for (int i = 0; i < nz; ++i)
        set(row[i], column[i], entry[i]);
    } finally {
      r.close();
    }
  }


  @Override
  public Matrix zero() {
    rows = new Linked();
    return this;
  }

  @Override
  public double get(int row, int column) {
    return rows.get(row, column);
  }

  @Override
  public void set(int row, int column, double value) {
    check(row, column);
    rows.set(row, column, value);
  }

  // avoids object creation
  static class ReusableMatrixEntry implements MatrixEntry {

    int row, col;
    double val;

    @Override
    public int column() {
      return col;
    }

    @Override
    public int row() {
      return row;
    }

    @Override
    public double get() {
      return val;
    }

    @Override
    public void set(double value) {
      throw new UnsupportedOperationException();
    }

    @Override
    public String toString() {
      return row + "," + col + "=" + val;
    }
  }

  @Override
  public Iterator<MatrixEntry> iterator() {
    return new Iterator<MatrixEntry>() {
      Node cur = rows.node;
      ReusableMatrixEntry entry = new ReusableMatrixEntry();

      @Override
      public boolean hasNext() {
        return cur != null;
      }

      @Override
      public MatrixEntry next() {
        entry.row = cur.row;
        entry.col = cur.col;
        entry.val = cur.val;
        cur = cur.tail;
        return entry;
      }

      @Override
      public void remove() {
        throw new UnsupportedOperationException("TODO");
      }
    };
  }

  @Override
  public Vector multAdd(double alpha, Vector x, Vector y) {
    return super.multAdd(alpha, x, y);
  }

  @Override
  public Vector transMultAdd(double alpha, Vector x, Vector y) {
    return super.transMultAdd(alpha, x, y);
  }

  @Override
  public Matrix multAdd(double alpha, Matrix B, Matrix C) {
    return super.multAdd(alpha, B, C);
  }

  @Override
  public Matrix transAmultAdd(double alpha, Matrix B, Matrix C) {
    return super.transAmultAdd(alpha, B, C);
  }

  @Override
  public Matrix transBmultAdd(double alpha, Matrix B, Matrix C) {
    return super.transBmultAdd(alpha, B, C);
  }

  @Override
  public Matrix transABmultAdd(double alpha, Matrix B, Matrix C) {
    return super.transABmultAdd(alpha, B, C);
  }
}
