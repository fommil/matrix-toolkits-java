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

    Node head;

    Node[] rows = new Node[numRows];

    // true if the head exists and has this row/col
    private boolean isHead(int row, int col) {
      return head != null && head.row == row && head.col == col;
    }

    // true if node exists, it's tail exists, and has this row/col
    private boolean isNext(Node node, int row, int col) {
      return node != null && node.tail != null && node.tail.row == row && node.tail.col == col;
    }

    public double get(int row, int col) {
      if (isHead(row, col))
        return head.val;
      Node node = findPreceeding(row, col);
      if (isNext(node, row, col))
        return node.tail.val;
      else
        return 0;
    }

    public void set(int row, int col, double val) {
      if (val == 0) {
        delete(row, col);
        return;
      }
      if (isHead(row, col))
        head.val = val;
      else {
        Node node = findPreceeding(row, col);
        if (node == null) {
          head = new Node(row, col, val, head);
          updateCache(head);
        } else if (isNext(node, row, col))
          node.tail.val = val;
        else {
          node.tail = new Node(row, col, val, node.tail);
          updateCache(node.tail);
        }
      }
    }

    private void updateCache(Node node) {
      if (rows[node.row] == null || node.col > rows[node.row].col)
        rows[node.row] = node;
    }

    private void delete(int row, int col) {
      Node node = findPreceeding(row, col);
      if (isHead(row, col)) {
        if (rows[row] == head) rows[row] = null;
        head = head.tail;
      } else if (isNext(node, row, col)) {
        if (rows[row] == node.tail)
          rows[row] = node.row == row ? node : null;
        node.tail = node.tail.tail;
      }
    }

    // returns the node that either references this
    // index, or should reference it if inserted.
    // null if there are no entries.
    Node findPreceeding(int row, int col) {
      Node last = row > 0 ? cached(row - 1) : null;
      Node cur = last != null ? last : head;
      while (cur != null && cur.row <= row) {
        if (cur.row == row && cur.col >= col) return last;
        last = cur;
        cur = cur.tail;
      }
      return last;
    }

    // helper for findPreceeding
    private Node cached(int row) {
      for (int i = row; i > 0; i--) {
        Node cached = rows[row - 1];
        if (cached != null)
          return cached;
      }
      return null;
    }

    Node row(int row) {
      Node prec = findPreceeding(row, 0);
      if (prec != null && prec.tail != null && prec.tail.row == row)
        return prec.tail;
      if (head != null && head.row == row)
        return head;
      return null;
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
      Node cur = rows.head;
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
  public Matrix scale(double alpha) {
    if (alpha == 0)
      zero();
    else if (alpha != 1) for (MatrixEntry e : this)
      set(e.row(), e.column(), e.get() * alpha);
    return this;
  }

  @Override
  public Matrix copy() {
    return new LinkedSparseMatrix(this);
  }

  @Override
  public Matrix transpose() {
    Linked old = rows;
    numRows = numColumns;
    numColumns = old.rows.length;
    rows = new Linked();
    Node node = old.head;
    while (node != null) {
      set(node.col, node.row, node.val);
      node = node.tail;
    }
    return this;
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
    checkMultAdd(B, C);
    for (int i = 0; i < numRows; i++) {
      Node row = rows.row(i);
      // TODO: optimise based on RHS Matrix
      if (row != null && alpha != 0)
        for (int j = 0; j < B.numColumns(); j++) {
          Node node = row;
          double v = 0;
          while (node != null && node.row == i) {
            v += (B.get(node.col, j) * node.val);
            node = node.tail;
          }
          if (v != 0) C.add(i, j, alpha * v);
        }
    }
    return C;
  }

  private double[] getCol(Matrix B, int col) {
    double[] array = new double[B.numRows()];
    for (int i = 0; i < B.numRows(); i++)
      array[i] = B.get(i, col);
    return array;
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
