package nl.bio.inf.peptidomicswebapp.models;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Plot Test class
 * @author Wouter Zeevat
 */
class PlotTest {

    static String path;

    @BeforeAll
    static void setup() {
        Path resourceDirectory = Paths.get("src", "test", "resources");
        path = resourceDirectory.toFile().getAbsolutePath();
    }

    @Test
    void createPlot() throws IOException {
        Plot plot = new Plot(Arrays.toString(Files.readAllBytes(Paths.get(path + "/1b58.pdb"))));
        assertNotNull(plot);
    }

    @Test
    void bytes() throws IOException {
        Plot plot = new Plot(Arrays.toString(Files.readAllBytes(Paths.get(path + "/1b58.pdb"))));
        assertEquals(Arrays.toString(Files.readAllBytes(Paths.get(path + "/1b58.pdb"))), plot.bytes());
    }

}