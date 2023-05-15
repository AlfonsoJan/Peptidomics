package nl.bio.inf.peptidomicswebapp.service;

/**
 * Interface for the pythonService class.
 * @author Wouter Zeevat
 */
public interface PythonConstructor {

    String getChainsPBD(String pythonPath, String pdbPath);

    String PDBAnalyse(String pythonPath, String pdbPath, String param, String pdbCode);
}