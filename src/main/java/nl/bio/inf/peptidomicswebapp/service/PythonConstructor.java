package nl.bio.inf.peptidomicswebapp.service;

public interface PythonConstructor {

    void createTempNumpyFile(String pythonPath, String uniqueNameNumpy, String pdbPath, String parameter);
    String createPcaPlot(String pythonPath, String numpyPath);
    String createScatterPlot(String pythonPath, String numpyPath);
    String createPlotlyPcaPlot(String pythonPath, String numpyPath);
    String getChainsPBD(String pythonPath, String pdbID);

    String PDBAnalyse(String pythonPath, String pdbID, String param, String comparePDB);
}