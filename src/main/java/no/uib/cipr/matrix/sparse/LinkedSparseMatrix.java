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
 * a good fit for unstructured sparse matrices. A secondary
 * link maintains fast transpose iteration.
 * <p/>
 * However, memory requirements
 * ({@code 1 instance (64 bits), 2 int (2 x 64 bits), 2 ref (2 x 64 bits), 1 double (128 bits) = 448 bits}
 * per matrix element, plus {@code 2.n ref}s for the cache) are slightly higher
 * than structured sparse matrix storage.
 *
 * @author Sam Halliday
 */
@Log
public class LinkedSparseMatrix extends AbstractMatrix {

  // java.util.LinkedList is doubly linked and therefore too heavyweight.
  @AllArgsConstructor
  @ToString(exclude = {"rowTail", "colTail"})
  static class Node {
    final int row, col;
    double val;
    Node rowTail, colTail;
  }

  // there is a lot of duplicated code in this class between
  // row and col linkages, but subtle differences make it
  // extremely difficult to factor away.
  class Linked {

    Node headRow, headCol;

    Node[] rows = new Node[numRows], cols = new Node[numColumns];

    // true if the head exists and has this row/col
    private boolean isHeadRow(int row, int col) {
      return headRow != null && headRow.row == row && headRow.col == col;
    }

    private boolean isHeadCol(int row, int col) {
      return headCol != null && headCol.col == col && headCol.row == row;
    }

    // true if node exists, it's tail exists, and has this row/col
    private boolean isNextByRow(Node node, int row, int col) {
      return node != null && node.rowTail != null && node.rowTail.row == row && node.rowTail.col == col;
    }

    private boolean isNextByCol(Node node, int row, int col) {
      return node != null && node.colTail != null && node.colTail.col == col && node.colTail.row == row;
    }

    public double get(int row, int col) {
      if (isHeadRow(row, col))
        return headRow.val;
      if (isHeadCol(row, col))
        return headCol.val;
      if (col <= row) {
        Node node = findPreceedingByRow(row, col);
        if (isNextByRow(node, row, col))
          return node.rowTail.val;
      } else {
        Node node = findPreceedingByCol(row, col);
        if (isNextByCol(node, row, col))
          return node.colTail.val;
      }
      return 0;
    }

    public void set(int row, int col, double val) {
      if (val == 0) {
        delete(row, col);
        return;
      }
      if (isHeadRow(row, col))
        headRow.val = val;
      else {
        Node node = findPreceedingByRow(row, col);
        if (node == null) {
          headRow = new Node(row, col, val, headRow, headRow != null ? headRow.colTail : null);
          updateRowCache(headRow);
        } else if (isNextByRow(node, row, col))
          node.rowTail.val = val;
        else {
          node.rowTail = new Node(row, col, val, node.rowTail, node.colTail);
          updateRowCache(node.rowTail);
        }
      }
      if (isHeadCol(row, col))
        headCol.val = val;
      else {
        Node node = findPreceedingByCol(row, col);
        if (node == null) {
          headCol = new Node(row, col, val, headCol != null ? headCol.rowTail : null, headCol);
          updateColCache(headCol);
        } else if (isNextByCol(node, row, col))
          node.colTail.val = val;
        else {
          node.colTail = new Node(row, col, val, node.rowTail, node.colTail);
          updateColCache(node.colTail);
        }
      }
    }

    private void updateRowCache(Node inserted) {
      if (rows[inserted.row] == null || inserted.col > rows[inserted.row].col)
        rows[inserted.row] = inserted;
    }

    private void updateColCache(Node inserted) {
      if (cols[inserted.col] == null || inserted.row > cols[inserted.col].row)
        cols[inserted.col] = inserted;
    }

    private void delete(int row, int col) {
      Node node = findPreceedingByRow(row, col);
      if (isHeadRow(row, col)) {
        if (rows[row] == headRow) rows[row] = null;
        headRow = headRow.rowTail;
      } else if (isNextByRow(node, row, col)) {
        if (rows[row] == node.rowTail)
          rows[row] = node.row == row ? node : null;
        node.rowTail = node.rowTail.rowTail;
      }
    }

    // returns the node that either references this
    // index, or should reference it if inserted.
    // null if there are no entries.
    Node findPreceedingByRow(int row, int col) {
      // FIXME: check the -1 here
      Node last = row > 0 ? cachedByRow(row - 1) : null;
      Node cur = last != null ? last : headRow;
      while (cur != null && cur.row <= row) {
        if (cur.row == row && cur.col >= col) return last;
        last = cur;
        cur = cur.rowTail;
      }
      return last;
    }

    // helper for findPreceeding
    private Node cachedByRow(int row) {
      for (int i = row; i > 0; i--) {
        Node cached = rows[row - 1];
        if (cached != null)
          return cached;
      }
      return null;
    }

    Node findPreceedingByCol(int row, int col) {
      // FIXME: check the -1 here
      Node last = col > 0 ? cachedByCol(col - 1) : null;
      Node cur = last != null ? last : headCol;
      while (cur != null && cur.col <= col) {
        if (cur.col == col && cur.row >= row) return last;
        last = cur;
        cur = cur.colTail;
      }
      return last;
    }

    private Node cachedByCol(int col) {
      for (int i = col; i > 0; i--) {
        Node cached = cols[col - 1];
        if (cached != null)
          return cached;
      }
      return null;
    }


    Node startOfRow(int row) {
      Node prec = findPreceedingByRow(row, 0);
      if (prec != null && prec.rowTail != null && prec.rowTail.row == row)
        return prec.rowTail;
      if (headRow != null && headRow.row == row)
        return headRow;
      return null;
    }

    Node startOfCol(int col) {
      Node prec = findPreceedingByCol(0, col);
      if (prec != null && prec.colTail != null && prec.colTail.col == col)
        return prec.colTail;
      if (headCol != null && headCol.col == col)
        return headCol;
      return null;
    }
  }

  Linked links;

  protected LinkedSparseMatrix(int numRows, int numColumns) {
    super(numRows, numColumns);
    links = new Linked();
  }

  protected LinkedSparseMatrix(Matrix A) {
    super(A);
    links = new Linked();
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
      links = new Linked();

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
    links = new Linked();
    return this;
  }

  @Override
  public double get(int row, int column) {
    return links.get(row, column);
  }

  @Override
  public void set(int row, int column, double value) {
    check(row, column);
    links.set(row, column, value);
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
      Node cur = links.headRow;
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
        cur = cur.rowTail;
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
    Linked old = links;
    numRows = numColumns;
    numColumns = old.rows.length;
    links = new Linked();
    Node node = old.headRow;
    while (node != null) {
      set(node.col, node.row, node.val);
      node = node.rowTail;
    }
    return this;
  }

  @Override
  public Vector multAdd(double alpha, Vector x, Vector y) {
    checkMultAdd(x, y);
    if (alpha == 0) return y;
    Node node = links.headRow;
    while (node != null) {
      y.add(node.row, alpha * node.val * x.get(node.col));
      node = node.rowTail;
    }
    return y;
  }

  @Override
  public Vector transMultAdd(double alpha, Vector x, Vector y) {
    checkTransMultAdd(x, y);
    if (alpha == 0) return y;
    Node node = links.headCol;
    while (node != null) {
      y.add(node.col, alpha * node.val * x.get(node.row));
      node = node.colTail;
    }
    return y;
  }

  // TODO: optimise matrix mults based on RHS Matrix

  @Override
  public Matrix multAdd(double alpha, Matrix B, Matrix C) {
    checkMultAdd(B, C);
    if (alpha == 0)
      return C;
    for (int i = 0; i < numRows; i++) {
      Node row = links.startOfRow(i);
      if (row != null)
        for (int j = 0; j < B.numColumns(); j++) {
          Node node = row;
          double v = 0;
          while (node != null && node.row == i) {
            v += (B.get(node.col, j) * node.val);
            node = node.rowTail;
          }
          if (v != 0) C.add(i, j, alpha * v);
        }
    }
    return C;
  }

  @Override
  public Matrix transBmultAdd(double alpha, Matrix B, Matrix C) {
    checkTransBmultAdd(B, C);
    if (alpha == 0)
      return C;
    for (int i = 0; i < numRows; i++) {
      Node row = links.startOfRow(i);
      if (row != null)
        for (int j = 0; j < B.numRows(); j++) {
          Node node = row;
          double v = 0;
          while (node != null && node.row == i) {
            v += (B.get(j, node.col) * node.val);
            node = node.rowTail;
          }
          if (v != 0) C.add(i, j, alpha * v);
        }
    }
    return C;
  }

  @Override
  public Matrix transAmultAdd(double alpha, Matrix B, Matrix C) {
    checkTransAmultAdd(B, C);
    if (alpha == 0)
      return C;
    for (int i = 0; i < numColumns; i++) {
      Node row = links.startOfCol(i);
      if (row != null)
        for (int j = 0; j < B.numColumns(); j++) {
          Node node = row;
          double v = 0;
          while (node != null && node.col == i) {
            v += (B.get(node.row, j) * node.val);
            node = node.colTail;
          }
          if (v != 0) C.add(i, j, alpha * v);
        }
    }
    return C;
  }

  @Override
  public Matrix transABmultAdd(double alpha, Matrix B, Matrix C) {
    checkTransABmultAdd(B, C);
    if (alpha == 0)
      return C;
    for (int i = 0; i < numColumns; i++) {
      Node row = links.startOfCol(i);
      if (row != null)
        for (int j = 0; j < B.numRows(); j++) {
          Node node = row;
          double v = 0;
          while (node != null && node.col == i) {
            v += (B.get(j, node.row) * node.val);
            node = node.colTail;
          }
          if (v != 0) C.add(i, j, alpha * v);
        }
    }
    return C;
  }
}
