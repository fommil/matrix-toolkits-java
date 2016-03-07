package no.uib.cipr.matrix.io;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixTestAbstract;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class MatrixVectorIoTest {

    @Test
    public void testWriteRead() throws Exception {
        DenseMatrix mat = new DenseMatrix(
                new double[][]{{1.1, 1.2}, {1.3, 1.4}});
        File matrixFile = new File("TestMatrixFile");
        BufferedWriter out = new BufferedWriter(new FileWriter(matrixFile));
        MatrixVectorWriter writer = new MatrixVectorWriter(out);
        MatrixInfo mi = new MatrixInfo(false, MatrixInfo.MatrixField.Real,
                MatrixInfo.MatrixSymmetry.General);
        writer.printMatrixInfo(mi);
        writer.printMatrixSize(new MatrixSize(mat.numRows(), mat.numColumns(),
                mat.numColumns() * mat.numRows()), mi);
        writer.printArray(mat.getData());
        writer.close();
        Matrix newMat = new DenseMatrix(new MatrixVectorReader(new FileReader(
                matrixFile)));

        MatrixTestAbstract.assertMatrixEquals(mat, newMat);
    }

    @Test
    public void testSparseWriteRead() throws Exception {
        CompRowMatrix mat = new CompRowMatrix(3, 3, new int[][]{{1, 2}, {0, 2},
                {1}});
        mat.set(0, 1, 1);
        mat.set(0, 2, 2);
        mat.set(1, 0, 3);
        mat.set(1, 2, 4);
        mat.set(2, 1, 5);
        File testFile = new File("TestMatrixFile");
        testFile.deleteOnExit();
        BufferedWriter out = new BufferedWriter(new FileWriter(testFile));
        MatrixVectorWriter writer = new MatrixVectorWriter(out);
        MatrixInfo mi = new MatrixInfo(true, MatrixInfo.MatrixField.Real,
                MatrixInfo.MatrixSymmetry.General);
        writer.printMatrixInfo(mi);
        writer.printMatrixSize(
                new MatrixSize(mat.numColumns(), mat.numColumns(), mat
                        .getData().length), mi);
        int[] rows = buildRowArray(mat);
        writer.printCoordinate(rows, mat.getColumnIndices(), mat.getData(), 1);
        out.close();
        CompRowMatrix newMat = new CompRowMatrix(new MatrixVectorReader(
                new FileReader(testFile)));
        MatrixTestAbstract.assertMatrixEquals(mat, newMat);
    }

    private int[] buildRowArray(CompRowMatrix mat) {
        int[] rows = new int[mat.getData().length];
        int curRow = -1;
        int nextValidRow = 0;
        int curIndex = 0;
        int[] rowPs = mat.getRowPointers();
        // find first valid row
        while (nextValidRow < mat.numRows() && rowPs[nextValidRow] < 0)
            nextValidRow++;
        while (true) {
            curRow = nextValidRow;
            nextValidRow++;
            while (nextValidRow < mat.numRows() && rowPs[nextValidRow] < 0)
                nextValidRow++;
            // enter the remainder of data as currentRow
            if (nextValidRow >= mat.numRows()) {
                while (curIndex < mat.getData().length)
                    rows[curIndex++] = curRow;
                break;
            }
            int nextRowIndex = rowPs[nextValidRow];
            while (curIndex < nextRowIndex)
                rows[curIndex++] = curRow;
        }
        return rows;
    }
}
