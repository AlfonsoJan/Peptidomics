package nl.bio.inf.peptidomicswebapp.service;

public interface PythonConstructor {

    String getChainsPBD(String pythonPath, String pdbID);

    String PDBAnalyse(String pythonPath, String pdbID, String param, String comparePDB);
}