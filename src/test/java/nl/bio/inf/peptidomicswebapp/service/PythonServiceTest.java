package nl.bio.inf.peptidomicswebapp.service;

import nl.bio.inf.peptidomicswebapp.exceptions.InvalidPDBCodeException;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.springframework.boot.test.context.SpringBootTest;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.*;

/**
 * PythonService Test class
 * @author Wouter Zeevat
 */
@SpringBootTest
class PythonServiceTest {

    private final PythonService pythonService = new PythonService();

    static String testResources;
    static String resources;

    @BeforeAll
    static void setup() {
        Path resourceDirectory = Paths.get("src", "test", "resources");
        testResources = resourceDirectory.toFile().getAbsolutePath();

        resourceDirectory = Paths.get("data");
        resources = resourceDirectory.toFile().getAbsolutePath();
    }

    @Test
    @DisplayName("Tests the get chains method with a normal file and pdb file")
    void testGetChainsPDB() {
        String pythonDir = resources + "/scripts/retrieve_chains_pdb.py";
        String PDBDir = testResources + "/1b58.pdb";

        String chain = pythonService.getChainsPBD(pythonDir, PDBDir);
        assert(chain.equalsIgnoreCase("{\"A\": \"518\", \"B\": \"4\"}"));
    }

    @Test
    @DisplayName("Tests the get chains method with a normal file but null PDB!")
    void testGetChainsPDBNullPDB() {
        String pythonDir = resources + "/scripts/retrieve_chains_pdb.py";
        assertThrows(NullPointerException.class, () -> pythonService.getChainsPBD(pythonDir, null));
    }

    @Test
    @DisplayName("Tests the get chains method with a normal file but null PDB!")
    void testGetChainsPDBNullScript() {
        String PDBDir = testResources + "/1b58.pdb";
        assertThrows(NullPointerException.class, () -> pythonService.getChainsPBD(null, PDBDir));
    }

    @Test
    @DisplayName("Tests the get chains method with a invalid script")
    void testGetChainsPDBInvalidScript() {
        String pythonDir = "nothing";
        String PDBDir = testResources + "/1b58.pdb";
        assertThrows(RuntimeException.class, () -> pythonService.getChainsPBD(pythonDir, PDBDir));

    }

    @Test
    @DisplayName("Tests the get chains method with a invalid script")
    void testGetChainsPDBInvalidPDB() {
        String pythonDir = resources + "/scripts/retrieve_chains_pdb.py";
        String PDBDir = "nothing";
        assertThrows(RuntimeException.class, () -> pythonService.getChainsPBD(pythonDir, PDBDir));

    }


}