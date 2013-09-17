package no.uib.cipr.matrix.sparse;

import lombok.AllArgsConstructor;
import no.uib.cipr.matrix.AbstractMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.io.MatrixVectorReader;

import java.io.IOException;

/**
 * A Sorted Linked List implementation of Matrix that
 * has O(ln(n)) insertion cost with an ordering optimised
 * for matrix or vector multiplication: a good fit for
 * unstructured sparse matrices.
 * <p/>
 * However, transpose
 * multiplications are costly, individual entry lookups are expensive
 * and memory requirements are higher than the other structured sparse
 * matrix types.
 * <p/>
 * TODO: straight or transpose mult optimisation should be a user option... potentially both! (at memory cost)
 *
 * @author Sam Halliday
 */
public class LinkedListSparseMatrix extends AbstractMatrix {

  // the next 100 lines could be done in one line of Scala... *sigh*
  @AllArgsConstructor
  private static class Entry implements Comparable<Entry> {
    private final int row, col;
    private double val;

    @Override
    public int compareTo(Entry o) {
      throw new UnsupportedOperationException("TODO");
    }
  }

  // java.util.LinkedList is doubly linked
  // and therefore too heavyweight.
  @AllArgsConstructor
  private static class Node {
    private final Entry head;
    private Node tail;
  }

  private static class RowLinked {

    private Node node;

    // returns the node that either references this
    // index, or would be the one to update should it
    // be inserted: null if there are no entries.
    private Node findPreceeding(int row, int col) {
      assert row >= 0 && col >= 0;
      Node last = null;
      Node cur = node;
      while (cur != null && cur.head.row <= row && cur.head.col < col) {
        last = cur;
        cur = cur.tail;
      }
      return last;
    }

    public double get(int row, int col) {
      Node prec = findPreceeding(row, col);
      if (prec == null) return 0;
      if (prec.tail == null) return 0;
      if (prec.tail.head.row != row || prec.tail.head.col != col) return 0;
      return prec.tail.head.val;
    }

    // also functions as insert
    public void update(int row, int col, double val) {
      if (val == 0) {
        remove(row, col);
        return;
      }
      Node prec = findPreceeding(row, col);
      if (prec == null) {
        node = new Node(new Entry(row, col, val), null);
      } else if (prec.tail == null) {
        prec.tail = new Node(new Entry(row, col, val), null);
      } else if (prec.tail.head.row != row || prec.tail.head.col != col) {
        prec.tail = new Node(new Entry(row, col, val), prec.tail.tail);
      } else {
        prec.tail.head.val = val;
      }
    }

    public void remove(int row, int col) {
      Node prec = findPreceeding(row, col);
      if (prec == null) return;
      if (prec.tail == null) return;
      if (prec.tail.head.row != row || prec.tail.head.col != col) return;
      prec.tail = prec.tail.tail;
    }

  }


  protected LinkedListSparseMatrix(int numRows, int numColumns) {
    super(numRows, numColumns);
  }

  protected LinkedListSparseMatrix(Matrix A) {
    super(A);

    // TODO
  }

  public LinkedListSparseMatrix(MatrixVectorReader r) throws IOException {
    super(0, 0);
    // TODO
  }

}
