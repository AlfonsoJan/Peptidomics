package nl.bio.inf.peptidomicswebapp.service;

import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.springframework.core.io.ClassPathResource;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Objects;

import static org.junit.jupiter.api.Assertions.*;

class PythonServiceTest {

    private final PythonService service = new PythonService();

    static String testResources;

    @BeforeAll
    static void setup() {
        Path resourceDirectory = Paths.get("src", "test", "resources");
        testResources = resourceDirectory.toFile().getAbsolutePath();
    }

    @Test
    void createTempNumpyFile() throws Exception {

        File folderScripts = new ClassPathResource("scripts").getFile();
        File script = null;

        try {
            for (File f : Objects.requireNonNull(folderScripts.listFiles())) {
                if ("read_pdb.py".equals(f.getName())) {
                    script = f;
                }
            }

            service.createTempNumpyFile(
                    script.getAbsolutePath(),
                    testResources + "/output/pdb",
                    testResources.toString() + "/1b58.pdb",
                    "3");
        } catch (Exception e) {
            throw new Exception();
        }

        assertTrue(FileUtils.contentEquals(new File(testResources.toString() + "/output/pdb.npy"), new File(testResources + "/example.npy")));
    }

    @Test
    void createPcaPlot() throws Exception {

        File folderScripts = new ClassPathResource("scripts").getFile();
        File script = null;

        try {
            for (File f : Objects.requireNonNull(folderScripts.listFiles())) {
                if ("pca_dim_plot.py".equals(f.getName())) {
                    script = f;
                }
            }

            String bytes = service.createPcaPlot(
                    script.getAbsolutePath(),
                    testResources + "/example.npy");

            assertTrue(bytes.startsWith("iVBORw0KGgoAAAANSUhEUgAAAg0AAAF3CAYAAAAmSXiuAAAAOXRFWHRTb2Z0d2Fy")); // Start of the bytes example plot

        } catch (Exception e) {
            throw new Exception();
        }
    }

    @Test
    void createScatterPlot() throws Exception {

        File folderScripts = new ClassPathResource("scripts").getFile();
        File script = null;

        try {
            for (File f : Objects.requireNonNull(folderScripts.listFiles())) {
                if ("axis_scatter_plot.py".equals(f.getName())) {
                    script = f;
                }
            }

            String bytes = service.createScatterPlot(
                    script.getAbsolutePath(),
                    testResources + "/example.npy");
            assertTrue(bytes.startsWith("iVBORw0KGgoAAAANSUhEUgAAA4QAAAF3CAYAAAD914WlAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZl")); // Start of the bytes example plot

        } catch (Exception e) {
            throw new Exception();
        }

    }

    @Test
    void createPlotlyPcaPlot() throws Exception {

        File folderScripts = new ClassPathResource("scripts").getFile();
        File script = null;

        try {
            for (File f : Objects.requireNonNull(folderScripts.listFiles())) {
                if ("plotly_pca_plot.py".equals(f.getName())) {
                    script = f;
                }
            }

            String bytes = service.createScatterPlot(
                    script.getAbsolutePath(),
                    testResources + "/example.npy");
            assertTrue(bytes.startsWith("{\"x\": [4.238776206970215, 3.4193129539489746, 3.3998899459838867, 0.07000348716974258,")); // Start of the bytes example plot

        } catch (Exception e) {
            throw new Exception();
        }
    }
}