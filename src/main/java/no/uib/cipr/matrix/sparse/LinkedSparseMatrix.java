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
 * ({@code 1 instance (8 bytes), 2 int (16 bytes), 2 ref (16 bytes), 1 double (8 bytes) = 48 bytes}
 * per matrix element, plus {@code 8 x numcol + 8 x numrow bytes}s for the cache) are slightly higher
 * than structured sparse matrix storage. Note that on 32 bit JVMs, or on 64 bit JVMs
 * with <a href="https://wikis.oracle.com/display/HotSpotInternals/CompressedOops">CompressedOops</a>
 * enabled, references and ints only cost 4 bytes each, bringing the cost to 28 bytes per element.
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

    final Node head = new Node(0, 0, 0, null, null);

    Node[] rows = new Node[numRows], cols = new Node[numColumns];

    private boolean isHead(int row, int col) {
      return head.row == row && head.col == col;
    }

    // true if node exists, it's row tail exists, and has this row/col
    private boolean isNextByRow(Node node, int row, int col) {
      return node != null && node.rowTail != null && node.rowTail.row == row && node.rowTail.col == col;
    }

    // true if node exists, it's col tail exists, and has this row/col
    private boolean isNextByCol(Node node, int row, int col) {
      return node != null && node.colTail != null && node.colTail.col == col && node.colTail.row == row;
    }

    public double get(int row, int col) {
      if (isHead(row, col))
        return head.val;
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
      if (isHead(row, col))
        head.val = val;
      else {
        Node prevRow = findPreceedingByRow(row, col);
        if (isNextByRow(prevRow, row, col))
          prevRow.rowTail.val = val;
        else {
          Node prevCol = findPreceedingByCol(row, col);
          Node nextCol = findNextByCol(row, col);
          prevRow.rowTail = new Node(row, col, val, prevRow.rowTail, nextCol);
          prevCol.colTail = prevRow.rowTail;
          updateCache(prevRow.rowTail);
        }
      }
      // DEBUGGING
      if (isHead(row, col))
        assert head.val == val;
      else {
        Node node = findPreceedingByCol(row, col);
        assert node != null;
        assert node.colTail.val == val;
      }
    }

    private Node findNextByCol(int row, int col) {
      Node cur = cachedByCol(col - 1);
      while (cur != null) {
        if (row < cur.row && col <= cur.col || col < cur.col) return cur;
        cur = cur.colTail;
      }
      return cur;
    }

    private void updateCache(Node inserted) {
      if (rows[inserted.row] == null || inserted.col > rows[inserted.row].col)
        rows[inserted.row] = inserted;
      if (cols[inserted.col] == null || inserted.row > cols[inserted.col].row)
        cols[inserted.col] = inserted;
    }

    private void delete(int row, int col) {
      if (isHead(row, col)) {
        head.val = 0;
        return;
      }
      Node precRow = findPreceedingByRow(row, col);
      Node precCol = findPreceedingByCol(row, col);
      if (isNextByRow(precRow, row, col)) {
        if (rows[row] == precRow.rowTail)
          rows[row] = precRow.row == row ? precRow : null;
        precRow.rowTail = precRow.rowTail.rowTail;
      }
      if (isNextByCol(precCol, row, col)) {
        if (cols[col] == precCol.colTail)
          cols[col] = precCol.col == col ? precCol : null;
        precCol.colTail = precCol.colTail.colTail;
      }
    }

    // returns the node that either references this
    // index, or should reference it if inserted.
    Node findPreceedingByRow(int row, int col) {
      Node last = cachedByRow(row - 1);
      Node cur = last;
      while (cur != null && cur.row <= row) {
        if (cur.row == row && cur.col >= col) return last;
        last = cur;
        cur = cur.rowTail;
      }
      return last;
    }

    // helper for findPreceeding
    private Node cachedByRow(int row) {
      for (int i = row; i >= 0; i--)
        if (rows[i] != null)
          return rows[i];
      return head;
    }

    Node findPreceedingByCol(int row, int col) {
      Node last = cachedByCol(col - 1);
      Node cur = last;
      while (cur != null && cur.col <= col) {
        if (cur.col == col && cur.row >= row) return last;
        last = cur;
        cur = cur.colTail;
      }
      return last;
    }

    private Node cachedByCol(int col) {
      for (int i = col; i >= 0; i--)
        if (cols[i] != null)
          return cols[i];
      return head;
    }

    Node startOfRow(int row) {
      if (row == 0) return head;
      Node prec = findPreceedingByRow(row, 0);
      if (prec.rowTail != null && prec.rowTail.row == row)
        return prec.rowTail;
      return null;
    }

    Node startOfCol(int col) {
      if (col == 0) return head;
      Node prec = findPreceedingByCol(0, col);
      if (prec != null && prec.colTail != null && prec.colTail.col == col)
        return prec.colTail;
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
      Node cur = links.head;
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
    Node node = old.head;
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
    Node node = links.head;
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
    Node node = links.head;
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
