/*
 * Copyright (C) 2003-2006 Bjørn-Ove Heimsund
 * 
 * This file is part of MTJ.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package no.uib.cipr.matrix.distributed;

/**
 * Performs a reduction operation. When performing a reduction, start with the
 * value return by the <i>init</i> function, for example:
 * 
 * <pre>
 * int[] x, y;
 * Reduction r;
 * // ...
 * r.initInt(x);
 * r.opInt(x, y);
 * </pre>
 * 
 * <p>
 * Many predefined reductions are available in
 * {@link no.uib.cipr.matrix.distributed.Reductions}.
 * </p>
 *
 * @deprecated the <code>no.uib.cipr.matrix.distributed</code> package has been deprecated because
 * of a number of hard to fix concurrency bugs. It is distributed only for backwards compatibility,
 * but is not recommended. The utility of this package is questionable, as it does not allow
 * distribution of computation between JVMs or across a network. For many people, distributed
 * computing of multiple matrices can be achieved at a user-level through the
 * <a href="http://jppf.org">JPPF Framework</a>.
 * Users who need to deal with few very large matrices may wish to implement their own storage classes
 * and solvers using JPPF, but this will not be supported directly in matrix-toolkits-java.
 */
@Deprecated
public abstract class Reduction {

    /**
     * Sets up the output data
     */
    public void init(Object x) {
        if (x instanceof double[])
            initDouble((double[]) x);
        else if (x instanceof int[])
            initInt((int[]) x);
        else if (x instanceof boolean[])
            initBoolean((boolean[]) x);
        else if (x instanceof byte[])
            initByte((byte[]) x);
        else if (x instanceof char[])
            initChar((char[]) x);
        else if (x instanceof short[])
            initShort((short[]) x);
        else if (x instanceof long[])
            initLong((long[]) x);
        else if (x instanceof float[])
            initFloat((float[]) x);
        else
            throw new IllegalArgumentException("Datatype is not supported");
    }

    /**
     * Adds to the output data
     * 
     * @param x
     *            Output data
     * @param y
     *            New input data
     */
    public void op(Object x, Object y) {
        if (x instanceof double[])
            opDouble((double[]) x, (double[]) y);
        else if (x instanceof int[])
            opInt((int[]) x, (int[]) y);
        else if (x instanceof boolean[])
            opBoolean((boolean[]) x, (boolean[]) y);
        else if (x instanceof byte[])
            opByte((byte[]) x, (byte[]) y);
        else if (x instanceof char[])
            opChar((char[]) x, (char[]) y);
        else if (x instanceof short[])
            opShort((short[]) x, (short[]) y);
        else if (x instanceof long[])
            opLong((long[]) x, (long[]) y);
        else if (x instanceof float[])
            opFloat((float[]) x, (float[]) y);
        else
            throw new IllegalArgumentException("Datatype is not supported");
    }

    protected abstract void initBoolean(boolean[] x);

    protected abstract void initByte(byte[] x);

    protected abstract void initChar(char[] x);

    protected abstract void initShort(short[] x);

    protected abstract void initInt(int[] x);

    protected abstract void initFloat(float[] x);

    protected abstract void initLong(long[] x);

    protected abstract void initDouble(double[] x);

    protected abstract void opBoolean(boolean[] x, boolean[] y);

    protected abstract void opByte(byte[] x, byte[] y);

    protected abstract void opChar(char[] x, char[] y);

    protected abstract void opShort(short[] x, short[] y);

    protected abstract void opInt(int[] x, int[] y);

    protected abstract void opFloat(float[] x, float[] y);

    protected abstract void opLong(long[] x, long[] y);

    protected abstract void opDouble(double[] x, double[] y);

}