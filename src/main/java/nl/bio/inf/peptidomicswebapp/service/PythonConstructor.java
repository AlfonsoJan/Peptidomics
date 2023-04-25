package nl.bio.inf.peptidomicswebapp.service;

/**
@Wouter Zeevat
 **/
public interface PythonConstructor {

    String getChainsPBD(String pythonPath, String pdbID);

    String PDBAnalyse(String pythonPath, String pdbID, String param, String comparePDB);
}