package nl.bio.inf.peptidomicswebapp.models;

import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.springframework.mock.web.MockMultipartFile;
import org.springframework.web.multipart.MultipartFile;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import static org.junit.jupiter.api.Assertions.*;

class PDBTest {

    static String path;

    @BeforeAll
    static void setup() {
        Path resourceDirectory = Paths.get("src","test","resources");
        path = resourceDirectory.toFile().getAbsolutePath();
    }

    @Test
    void createPDBFromCode() throws IOException {
        PDB pdb = new PDB("6zdh");
        assertEquals(pdb.getFileName(), "6zdh.pdb");
    }

    @Test
    void createPDBFromFile() throws IOException {
        Path p = Paths.get(path + "/6zdh.pdb");
        String name = "6zdh.pdb";
        String originalFileName = "6zdh.pdb";
        String contentType = "text/plain";
        byte[] content = null;
        try {
            content = Files.readAllBytes(p);
        } catch (final IOException e) {
        }
        MultipartFile result = new MockMultipartFile(name,
                originalFileName, contentType, content);

        PDB pdb = new PDB(name, result.getBytes(), PDB.getStructureFromInputstream(result.getInputStream()));
        assertEquals(pdb.getFileName(), "6zdh.pdb");
    }

    @Test
    void getStructureFromInputstream() {
    }

    @Test
    void getBytes() throws IOException {
        PDB pdb = new PDB("1b58");

        Path p = Paths.get(path + "/1b58_generated.pdb");
        Files.write(p, pdb.getBytes());
        File generated = new File(p.toString());
        File old = new File(path + "/1b58.pdb");

        assertTrue(FileUtils.contentEquals(generated, old));
    }

    @ParameterizedTest
    @ValueSource(strings = {"6zdh", "1b58"})
    void getFileName(String code) throws IOException {
        PDB pdb = new PDB(code);
        assertEquals(pdb.getFileName(), code + ".pdb");
    }

    @ParameterizedTest
    @ValueSource(strings = {"6zdh", "1b58"})
    void getStructureId(String code) throws IOException {
        PDB pdb = new PDB(code);
        assertEquals(code, pdb.getStructureId());
    }
}