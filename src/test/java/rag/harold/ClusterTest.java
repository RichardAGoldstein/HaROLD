package rag.harold;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

public class ClusterTest {
    private Cluster cluster = null;
    private File resourcesDir = null;
    private GammaCalc gamma = GammaCalc.get(0);

    private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
    private final ByteArrayOutputStream errContent = new ByteArrayOutputStream();

    private static final String TEST_RESOURCES_PATH = "src/test/resources";
    private static final String DATA_DIRECTORY = "SmallData";
    private static final String DATA_FILE_LIST_NAME = "FileList_AAA";

    @Before
    public void setUp() {
        this.resourcesDir = new File(TEST_RESOURCES_PATH);
        System.setOut(new PrintStream(outContent));
        System.setErr(new PrintStream(errContent));
    }

    @After
    public void restoreStreams() {
        System.setOut(System.out);
        System.setErr(System.err);
    }

    private File getDataFiles() throws Exception {
        String fileList = resourcesDir.getAbsolutePath() + "/" + DATA_DIRECTORY + "/" + DATA_FILE_LIST_NAME;
        List<String> dataFiles = Files.readAllLines(Paths.get(fileList));

        File tmpFile = File.createTempFile("prefix", null);
        BufferedWriter tmpWriter = new BufferedWriter(new FileWriter(tmpFile));
        for (String line : dataFiles) {
            tmpWriter.write(resourcesDir + "/" + DATA_DIRECTORY + "/" + line + "\n");
        }
        tmpWriter.close();

        return tmpFile;
    }

    @Test
    public void clusterInstantiate() throws Exception {
        Options o = new Options();
        File dataFile = getDataFiles();
        long seed = 1234L;
        cluster = new Cluster(dataFile, 3, o.initialAlphaParams, this.gamma, seed, false);
        boolean deleted = dataFile.delete();
    }

    @Test
    public void clusterRun() throws Exception {
        Options o = new Options();
        File dataFile = getDataFiles();
        long seed = 1234L;
        cluster = new Cluster(dataFile, 3, o.initialAlphaParams, this.gamma, seed, false);
        cluster.run();
        boolean deleted = dataFile.delete();
    }
}
