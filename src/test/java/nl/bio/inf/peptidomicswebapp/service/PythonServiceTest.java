package nl.bio.inf.peptidomicswebapp.service;

import nl.bio.inf.peptidomicswebapp.exceptions.InvalidPDBCodeException;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.*;

/**
 * PythonService Test class
 * @author Wouter Zeevat
 */
class PythonServiceTest {

    private final PythonService pythonService = new PythonService();

    static String testResources;
    static String resources;

    @BeforeAll
    static void setup() {
        Path resourceDirectory = Paths.get("src", "test", "resources");
        testResources = resourceDirectory.toFile().getAbsolutePath();

        resourceDirectory = Paths.get("src", "main", "resources");
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

    @ParameterizedTest
    @ValueSource(ints = {0, 1, 2, 3})
    @DisplayName("Tests the analysis with one null each time")
    void testAnalyseIsNull(int index) {
        String pythonDir = resources + "/scripts/pdb_analysis.py";
        String PDBDir = testResources + "/1b58.pdb";
        ArrayList<String> params = new ArrayList<>();
        params.add(pythonDir);
        params.add(PDBDir);
        params.add("3");
        params.add("1b58");
        params.set(index, null);
        assertThrows(NullPointerException.class, () -> pythonService.PDBAnalyse(params.get(0), params.get(1), params.get(2), params.get(3)));
    }

    @ParameterizedTest
    @ValueSource(strings = {"none", "0", "31"})
    @DisplayName("Tests the analysis with an invalid pepsizes")
    void testAnalyseInvalidPepSize(String size) {
        String pythonDir = resources + "/scripts/pdb_analysis.py";
        String PDBDir = testResources + "/1b58.pdb";
        assertThrows(NumberFormatException.class, () -> pythonService.PDBAnalyse(pythonDir, PDBDir, size, "1b58"));
    }

    @ParameterizedTest
    @ValueSource(strings = {"131", "0", "31f1"})
    @DisplayName("Tests the analysis with an invalid pepsizes")
    void testAnalyseInvalidCode(String code) {
        String pythonDir = resources + "/scripts/pdb_analysis.py";
        String PDBDir = testResources + "/1b58.pdb";
        assertThrows(InvalidPDBCodeException.class, () -> pythonService.PDBAnalyse(pythonDir, PDBDir, "3", code));
    }


}