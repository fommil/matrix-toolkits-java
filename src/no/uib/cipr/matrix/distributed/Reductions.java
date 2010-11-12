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
 * Contains predefined reductions
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
public class Reductions {

    private Reductions() {
        // No need to instantiate
    }

    public static Reduction sum() {
        return new Sum();
    }

    public static Reduction product() {
        return new Product();
    }

    public static Reduction max() {
        return new Max();
    }

    public static Reduction min() {
        return new Min();
    }

    public static Reduction and() {
        return new And();
    }

    public static Reduction or() {
        return new Or();
    }

    /**
     * Does a sum
     */
    private static class Sum extends NumericalReduction {

        @Override
        protected void initByte(byte[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 0;
        }

        @Override
        protected void initChar(char[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 0;
        }

        @Override
        protected void initShort(short[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 0;
        }

        @Override
        protected void initInt(int[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 0;
        }

        @Override
        protected void initFloat(float[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 0;
        }

        @Override
        protected void initLong(long[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 0;
        }

        @Override
        protected void initDouble(double[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 0;
        }

        @Override
        protected void opByte(byte[] x, byte[] y) {
            assert x.length == y.length;
            for (int i = 0; i < x.length; ++i)
                x[i] += y[i];
        }

        @Override
        protected void opChar(char[] x, char[] y) {
            assert x.length == y.length;
            for (int i = 0; i < x.length; ++i)
                x[i] += y[i];
        }

        @Override
        protected void opShort(short[] x, short[] y) {
            assert x.length == y.length;
            for (int i = 0; i < x.length; ++i)
                x[i] += y[i];
        }

        @Override
        protected void opInt(int[] x, int[] y) {
            assert x.length == y.length;
            for (int i = 0; i < x.length; ++i)
                x[i] += y[i];
        }

        @Override
        protected void opFloat(float[] x, float[] y) {
            assert x.length == y.length;
            for (int i = 0; i < x.length; ++i)
                x[i] += y[i];
        }

        @Override
        protected void opLong(long[] x, long[] y) {
            assert x.length == y.length;
            for (int i = 0; i < x.length; ++i)
                x[i] += y[i];
        }

        @Override
        protected void opDouble(double[] x, double[] y) {
            assert x.length == y.length;
            for (int i = 0; i < x.length; ++i)
                x[i] += y[i];
        }
    }

    /**
     * Multiplies the elements
     */
    private static class Product extends NumericalReduction {

        @Override
        protected void initByte(byte[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 1;
        }

        @Override
        protected void initChar(char[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 1;
        }

        @Override
        protected void initShort(short[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 1;
        }

        @Override
        protected void initInt(int[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 1;
        }

        @Override
        protected void initFloat(float[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 1;
        }

        @Override
        protected void initLong(long[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 1;
        }

        @Override
        protected void initDouble(double[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = 1;
        }

        @Override
        protected void opByte(byte[] x, byte[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] *= y[i];
        }

        @Override
        protected void opChar(char[] x, char[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] *= y[i];
        }

        @Override
        protected void opShort(short[] x, short[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] *= y[i];
        }

        @Override
        protected void opInt(int[] x, int[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] *= y[i];
        }

        @Override
        protected void opFloat(float[] x, float[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] *= y[i];
        }

        @Override
        protected void opLong(long[] x, long[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] *= y[i];
        }

        @Override
        protected void opDouble(double[] x, double[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] *= y[i];
        }
    }

    /**
     * Finds max
     */
    private static class Max extends NumericalReduction {

        @Override
        protected void initByte(byte[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Byte.MIN_VALUE;
        }

        @Override
        protected void initChar(char[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Character.MIN_VALUE;
        }

        @Override
        protected void initShort(short[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Short.MIN_VALUE;
        }

        @Override
        protected void initInt(int[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Integer.MIN_VALUE;
        }

        @Override
        protected void initFloat(float[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Float.MIN_VALUE;
        }

        @Override
        protected void initLong(long[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Long.MIN_VALUE;
        }

        @Override
        protected void initDouble(double[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Double.MIN_VALUE;
        }

        @Override
        protected void opByte(byte[] x, byte[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = (byte) Math.max(x[i], y[i]);
        }

        @Override
        protected void opChar(char[] x, char[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = (char) Math.max(x[i], y[i]);
        }

        @Override
        protected void opShort(short[] x, short[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = (short) Math.max(x[i], y[i]);
        }

        @Override
        protected void opInt(int[] x, int[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Math.max(x[i], y[i]);
        }

        @Override
        protected void opFloat(float[] x, float[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Math.max(x[i], y[i]);
        }

        @Override
        protected void opLong(long[] x, long[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Math.max(x[i], y[i]);
        }

        @Override
        protected void opDouble(double[] x, double[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Math.max(x[i], y[i]);
        }
    }

    /**
     * Finds min
     */
    private static class Min extends NumericalReduction {

        @Override
        protected void initByte(byte[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Byte.MAX_VALUE;
        }

        @Override
        protected void initChar(char[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Character.MAX_VALUE;
        }

        @Override
        protected void initShort(short[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Short.MAX_VALUE;
        }

        @Override
        protected void initInt(int[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Integer.MAX_VALUE;
        }

        @Override
        protected void initFloat(float[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Float.MAX_VALUE;
        }

        @Override
        protected void initLong(long[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Long.MAX_VALUE;
        }

        @Override
        protected void initDouble(double[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Double.MAX_VALUE;
        }

        @Override
        protected void opByte(byte[] x, byte[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = (byte) Math.min(x[i], y[i]);
        }

        @Override
        protected void opChar(char[] x, char[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = (char) Math.min(x[i], y[i]);
        }

        @Override
        protected void opShort(short[] x, short[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = (short) Math.min(x[i], y[i]);
        }

        @Override
        protected void opInt(int[] x, int[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Math.min(x[i], y[i]);
        }

        @Override
        protected void opFloat(float[] x, float[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Math.min(x[i], y[i]);
        }

        @Override
        protected void opLong(long[] x, long[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Math.min(x[i], y[i]);
        }

        @Override
        protected void opDouble(double[] x, double[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] = Math.min(x[i], y[i]);
        }
    }

    private static abstract class NumericalReduction extends Reduction {

        @Override
        protected void initBoolean(boolean[] x) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void opBoolean(boolean[] x, boolean[] y) {
            throw new UnsupportedOperationException();
        }
    }

    private static abstract class BooleanReduction extends Reduction {

        @Override
        protected void initByte(byte[] x) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void initChar(char[] x) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void initShort(short[] x) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void initInt(int[] x) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void initFloat(float[] x) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void initLong(long[] x) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void initDouble(double[] x) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void opByte(byte[] x, byte[] y) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void opChar(char[] x, char[] y) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void opShort(short[] x, short[] y) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void opInt(int[] x, int[] y) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void opFloat(float[] x, float[] y) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void opLong(long[] x, long[] y) {
            throw new UnsupportedOperationException();
        }

        @Override
        protected void opDouble(double[] x, double[] y) {
            throw new UnsupportedOperationException();
        }
    }

    /**
     * And boolean reduction
     */
    private static class And extends BooleanReduction {

        @Override
        protected void initBoolean(boolean[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = true;
        }

        @Override
        protected void opBoolean(boolean[] x, boolean[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] &= y[i];
        }
    }

    /**
     * Or boolean reduction
     */
    private static class Or extends BooleanReduction {

        @Override
        protected void initBoolean(boolean[] x) {
            for (int i = 0; i < x.length; ++i)
                x[i] = false;
        }

        @Override
        protected void opBoolean(boolean[] x, boolean[] y) {
            for (int i = 0; i < x.length; ++i)
                x[i] |= y[i];
        }
    }
}
