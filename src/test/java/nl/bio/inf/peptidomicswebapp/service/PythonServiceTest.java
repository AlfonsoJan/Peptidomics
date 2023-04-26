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

/**
 * PythonService Test class
 * @author Wouter Zeevat
 */
class PythonServiceTest {

    private final PythonService service = new PythonService();

    static String testResources;

    @BeforeAll
    static void setup() {
        Path resourceDirectory = Paths.get("src", "test", "resources");
        testResources = resourceDirectory.toFile().getAbsolutePath();
    }

}