package nl.bio.inf.peptidomicswebapp.models;

import nl.bio.inf.peptidomicswebapp.exceptions.InvalidPDBCodeException;
import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.*;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.springframework.mock.web.MockMultipartFile;
import org.springframework.web.multipart.MultipartFile;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * PDB Test class
 * @author Wouter Zeevat
 */
class PDBTest {

    static String path;
    static List<String> tempPaths = new ArrayList<>();

    @BeforeAll
    @DisplayName("Setting the source directory!")
    static void setup() {
        Path resourceDirectory = Paths.get("src","test","resources");
        path = resourceDirectory.toFile().getAbsolutePath();
    }

    @AfterAll
    @DisplayName("Deletes the temporary files created in the tests!")
    static void deleteTemp() throws IOException {
        for (String tempPath : tempPaths) {
            Files.deleteIfExists(Path.of(tempPath));
        }
    }

    @Test
    @DisplayName("Creates a pdb object from a code!")
    void createPDBFromCode() throws IOException, InvalidPDBCodeException {
        PDB pdb = new PDB("6zdh");
        assertEquals(pdb.getFileName(), "6zdh.pdb");
    }

    @Test
    @DisplayName("Creates a pdb object from a file!")
    void createPDBFromFile() throws IOException {
        Path p = Paths.get(path + "/6zdh.pdb");
        String name = "6zdh.pdb";
        String originalFileName = "6zdh.pdb";
        String contentType = "text/plain";
        byte[] content = null;
        try {
            content = Files.readAllBytes(p);
        } catch (final IOException ignored) {
        }
        MultipartFile result = new MockMultipartFile(name,
                originalFileName, contentType, content);

        PDB pdb = new PDB(result.getBytes(), name);
        assertEquals(pdb.getFileName(), "6zdh.pdb");
    }

    @Test
    @DisplayName("Creates a pdb object from a null object!")
    void testPDBWithNull() {
        try {
            new PDB(null);
        } catch (Exception exception) {
            assertEquals(exception.getClass(), NullPointerException.class);
        }
    }

    @Test
    @DisplayName("Creates a pdb object from an empty string!")
    void testPDBWithEmptyString() {
        assertThrows(InvalidPDBCodeException.class, () -> new PDB(""));
    }

    @Test
    @DisplayName("Creates a pdb object from a code that's too long!")
    void testPDBWithLongerCode() {
        assertThrows(InvalidPDBCodeException.class, () -> new PDB("ABCDE"));
    }

    @Test
    @DisplayName("Creates a pdb object from a code that's too short!")
    void testPDBWithInvalidCode() {
        assertThrows(InvalidPDBCodeException.class, () -> new PDB("ABCD"));
    }

    @Test
    @DisplayName("Testing the pdb with a shorter code than allowed!")
    void testPDBWithShorterCode() {
        assertThrows(InvalidPDBCodeException.class, () -> new PDB("ABC"));
    }

    @Disabled
    @Test
    @DisplayName("Tests the getBytes function!")
    void getBytes() throws IOException, InvalidPDBCodeException {
        PDB pdb = new PDB("1b58");

        Path p = Paths.get(path + "/output/1b58_generated.pdb");
        Files.write(p, pdb.getBytes());
        File generated = new File(p.toString());
        File old = new File(path + "/1b58.pdb");

        assertTrue(FileUtils.contentEquals(generated, old));
    }

    @ParameterizedTest
    @ValueSource(strings = {"6zdh", "1b58"})
    @DisplayName("Tests the get file name function!")
    void getFileName(String code) throws IOException, InvalidPDBCodeException {
        PDB pdb = new PDB(code);
        assertEquals(pdb.getFileName(), code + ".pdb");
    }

    @ParameterizedTest
    @ValueSource(strings = {"6zdh", "1b58"})
    @DisplayName("Test the get id function!")
    void getStructureId(String code) throws IOException, InvalidPDBCodeException {
        PDB pdb = new PDB(code);
        assertEquals(code, pdb.getStructureId());
    }

    @ParameterizedTest
    @ValueSource(strings = {"6zdh", "1b58"})
    @DisplayName("Tests the get structure function separately!")
    void testGetStructureFromInputstream(String code) throws IOException, InvalidPDBCodeException {
        PDB pdb = new PDB(code);
        String resultCode = PDB.getStructureFromInputStream(pdb.getBytes());
        assert(resultCode.equalsIgnoreCase(code));
    }

    @ParameterizedTest
    @ValueSource(strings = {"6zdh", "1b58"})
    @DisplayName("Tests the create temp file separately!")
    void testCreateTempFile(String code) throws IOException, InvalidPDBCodeException {
        PDB pdb = new PDB(code);
        String path = pdb.createTempFile();
        tempPaths.add(path);
        File file = new File(path);
        assertTrue(file.exists());
    }

    @Test
    @DisplayName("Tests the is valid on a normal pdb file from a code")
    void testIsValidCode() throws InvalidPDBCodeException, IOException {
        PDB pdb = new PDB("1b58");
        assertTrue(pdb.isValid());
    }

    @Test
    @DisplayName("Tests the is valid on a normal pdb file from a file")
    void testIsValidFile() throws IOException {
        Path p = Paths.get(path + "/6zdh.pdb");
        String name = "6zdh.pdb";
        String originalFileName = "6zdh.pdb";
        String contentType = "text/plain";
        byte[] content = null;
        try {
            content = Files.readAllBytes(p);
        } catch (final IOException ignored) {
        }
        MultipartFile result = new MockMultipartFile(name,
                originalFileName, contentType, content);

        PDB pdb = new PDB(result.getBytes(), name);
        assertTrue(pdb.isValid());
    }

    @Test
    @DisplayName("Tests the is valid on an invalid pdb file from a file")
    void testIsValidFileIsInvalid() throws IOException {
        Path p = Paths.get(path + "/invalid.pdb");
        String name = "invalid.pdb";
        String originalFileName = "invalid.pdb";
        String contentType = "text/plain";
        byte[] content = null;
        try {
            content = Files.readAllBytes(p);
        } catch (final IOException ignored) {
        }
        MultipartFile result = new MockMultipartFile(name,
                originalFileName, contentType, content);

        PDB pdb = new PDB(result.getBytes(), name);
        assertFalse(pdb.isValid());
    }
}