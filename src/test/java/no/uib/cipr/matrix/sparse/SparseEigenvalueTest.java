package no.uib.cipr.matrix.sparse;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.Vector;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import org.apache.commons.math3.linear.EigenDecomposition;

import java.util.Arrays;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;

public class SparseEigenvalueTest  {

    private CompColMatrix A;

    Random random;

    private static final int MIN_MATRIX_DIMENSION = 100;
    private static final int MAX_MATRIX_DIMENSION = 1000;

    /**
     * Creates a random sparse Matrix for eigenvalue computation
     * @throws Exception
     */
    @Before
    public void setUp() throws Exception {
        random = new Random();

        A = createRandomMatrix(MIN_MATRIX_DIMENSION, MAX_MATRIX_DIMENSION, 5);
    }

    private CompColMatrix createRandomMatrix(int minDim, int maxDim, int averageDensity) {
        if (averageDensity < 1){
            throw new IllegalArgumentException("Average density of Matrix must be >= 1");
        }
        int dim = minDim < maxDim ?
                minDim + random.nextInt(maxDim-minDim) :
                minDim;

        int[][] nz = new int[dim][];
        for (int col = 0; col < dim; col++){
            int valsPerCol = 1 + Math.min(dim-1, random.nextInt(averageDensity*2 - 1));
            nz[col] = new int[valsPerCol];

            TreeSet<Integer> targetRows = new TreeSet<Integer>();

            int fixedTarget = col+1>=dim?0:col+1;
            targetRows.add(fixedTarget); // create a circle

            for (int r = 0; r < valsPerCol-1; r++){
                int randomTarget = fixedTarget;
                while (targetRows.contains(randomTarget)){
                    randomTarget = random.nextInt(dim);
                }
                targetRows.add(randomTarget);
            }
            int r = 0;
            for (Integer target : targetRows){
                nz[col][r++] = target;
            }
        }

        CompColMatrix A = new CompColMatrix(dim, dim, nz);

        for (int c = 0; c < nz.length; c++){
            for (int r = 0; r < nz[c].length; r++){
                int row = nz[c][r];
                A.set(row,c,1+random.nextInt(5));
            }
        }
        return A;
    }

    @Test
    public void testEigenValue() throws Exception {
        ArpackGen generalSolver = new ArpackGen(A);
        generalSolver.setComputeOnlyEigenvalues(true);
        Map<Double, DenseVectorSub> eigenValueMap = generalSolver.solve(3, ArpackGen.Ritz.LR);

        RealMatrix commonsA = new Array2DRowRealMatrix(A.numRows(), A.numColumns());
        double[] data = A.getData();
        int col = 0;
        int row;
        double val;
        for (int elem = 0; elem < data.length; elem++){
            row = A.getRowIndices()[elem];
            if (A.getColumnPointers()[col+1] <= elem){
                col++;
            }
            val = A.get(row, col);
            commonsA.setEntry(row, col, val);

        }
        EigenDecomposition decomposition = new EigenDecomposition(commonsA);
        double[] realEigenvalues = decomposition.getRealEigenvalues();
        Arrays.sort(realEigenvalues);
        double largestEigenValue = realEigenvalues[realEigenvalues.length-1];

        double largestARPACKEigenValue = eigenValueMap.keySet().iterator().next();

        System.out.println("Apache: " + largestEigenValue+" vs. ARPACK: "+largestARPACKEigenValue);
        Assert.assertEquals(largestEigenValue,largestARPACKEigenValue, 10e-10);
    }


    /**
     * TODO: Gather EigenVectors correctly from ARPACK
     * @throws Exception
     */
    @Test
    @Ignore
    public void testEigenVectors() throws Exception {
        testMatrix(A);
    }

    private void testMatrix(CompColMatrix matrix) {
        ArpackGen generalSolver = new ArpackGen(matrix);
        generalSolver.setComputeOnlyEigenvalues(false);
        Map<Double, DenseVectorSub> eigenValueMap = generalSolver.solve(Math.min(matrix.numColumns()-2,5), ArpackGen.Ritz.LR);


        // test if the relationship A*x = lambda*x is true:
        if (eigenValueMap.isEmpty()){
            throw new IllegalStateException("No non-complex Eigenvalues found!");
        }
        Double largestEigenvalue = eigenValueMap.keySet().iterator().next();
        DenseVector v = eigenValueMap.get(largestEigenvalue).copy();
        DenseVector x = new DenseVector(v.size());
        matrix.mult(v, x);
        DenseVector vScaled = v.copy().scale(largestEigenvalue); // v is otherwise unit norm
        for (int i = 0; i < x.size(); i++){
            Assert.assertEquals(vScaled.get(i), x.get(i), 1e-7);
        }

        System.out.println("Eigenvalue "+largestEigenvalue+" matches with eigenvector:");

        System.out.println(v.toString());
    }

    @Test
    public void testExhaustively() throws Exception {
        int dim = 5;
        int errorcount = 0;
        int successCount = 0;
        int iteration = 0;
        DescriptiveStatistics errors = new DescriptiveStatistics();
        for (iteration = 0; iteration < 1000; iteration++){
            CompColMatrix m = createRandomMatrix(dim+iteration,dim+iteration,5);
            try {
                testMatrix(m);
                successCount++;
            } catch (IllegalStateException e){
                e.printStackTrace();
                errorcount++;
            } catch (AssertionError e) {
                String[] parts = e.getMessage().split("<");
                double expected = Double.valueOf(parts[1].split(">")[0]);
                double actual = Double.valueOf(parts[2].split(">")[0]);
                errors.addValue(expected-actual);
                errorcount++;
            }
        }
        System.out.println(successCount + " correct eigenvalue/eigenvector computations.");
        System.out.println(errorcount + " failures in "+ iteration+" iterations!");
        System.out.println("Convergence / imprecision errors "+errors.getN()+" \t average: "+errors.getMean());
    }
}
