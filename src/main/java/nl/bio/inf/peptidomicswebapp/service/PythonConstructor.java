package nl.bio.inf.peptidomicswebapp.service;

/**
 * Interface for the pythonService class.
 * @author Wouter Zeevat
 */
public interface PythonConstructor {

    String getChainsPBD(String pythonPath, String pdbID);

    String PDBAnalyse(String pythonPath, String pdbID, String param, String comparePDB);
}