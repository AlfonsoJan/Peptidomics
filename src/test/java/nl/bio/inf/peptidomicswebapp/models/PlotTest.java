package nl.bio.inf.peptidomicswebapp.models;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

class PlotTest {

    static String path;

    @BeforeAll
    static void setup() {
        Path resourceDirectory = Paths.get("src", "test", "resources");
        path = resourceDirectory.toFile().getAbsolutePath();
    }

    @Test
    void createPlot() throws IOException {
        Plot plot = new Plot(Arrays.toString(Files.readAllBytes(Paths.get(path + "/scatter_example.npy"))));
        assertNotNull(plot);
    }

    @Test
    void getBytes() throws IOException {
        Plot plot = new Plot(Arrays.toString(Files.readAllBytes(Paths.get(path + "/scatter_example.npy"))));
        assertEquals(Arrays.toString(Files.readAllBytes(Paths.get(path + "/scatter_example.npy"))), plot.getBytes());
    }

    @Test
    void setBytes() throws IOException {
        Plot plot = new Plot(Arrays.toString(Files.readAllBytes(Paths.get(path + "/scatter_example.npy"))));
        plot.setBytes("x");
        assertEquals("x", plot.getBytes());
    }
}