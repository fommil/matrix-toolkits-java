package no.uib.cipr.matrix.io;

import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixTestAbstract;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class MatrixVectorIoTest extends TestCase {

  public void testWriteRead() throws Exception {
    DenseMatrix mat = new DenseMatrix(new double[][] {{1.1,1.2},{1.3,1.4}});
    File matrixFile = new File("TestMatrixFile");
    BufferedWriter out = new BufferedWriter(new FileWriter(matrixFile));
    MatrixVectorWriter writer = new MatrixVectorWriter(out);
    writer.printMatrixInfo(new MatrixInfo(false, MatrixInfo.MatrixField.Real, MatrixInfo.MatrixSymmetry.General));
    writer.printMatrixSize(new MatrixSize(mat.numRows(), mat.numColumns(), mat.numColumns() * mat.numRows()));
    writer.printArray(mat.getData());
    writer.close();
    Matrix newMat = new DenseMatrix(new MatrixVectorReader(new FileReader(matrixFile)));

    MatrixTestAbstract.assertEquals(mat, newMat);
  }
}
