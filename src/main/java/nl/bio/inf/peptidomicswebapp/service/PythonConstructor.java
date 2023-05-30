package nl.bio.inf.peptidomicswebapp.service;

import nl.bio.inf.peptidomicswebapp.exceptions.InvalidPDBCodeException;

/**
 * Interface for the pythonService class.
 * @author Wouter Zeevat
 */
public interface PythonConstructor {

    String getChainsPBD(String pythonPath, String pdbPath);

    void PDBAnalyse(String pythonPath, String pdbPath, String param, String pdbCode) throws InvalidPDBCodeException;
}