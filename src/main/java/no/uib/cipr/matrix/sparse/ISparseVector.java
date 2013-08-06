/**
 * May 21, 2007
 * @author Samuel Halliday, ThinkTank Maths Limited
 * Copyright ThinkTank Maths Limited 2007
 */
package no.uib.cipr.matrix.sparse;

import no.uib.cipr.matrix.Vector;

/**
 * @author Samuel Halliday, ThinkTank Maths Limited
 */
public interface ISparseVector extends Vector {

	/**
	 * Returns the indices
	 */
	public int[] getIndex();

	/**
	 * Number of entries used in the sparse structure
	 */
	public int getUsed();
}
