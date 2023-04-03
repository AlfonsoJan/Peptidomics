package nl.bio.inf.peptidomicswebapp.service;

import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import org.springframework.stereotype.Service;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.logging.Logger;

@Service
public class PythonService implements PythonConstructor{
    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());
    private static final int EXIT_CODE = 0;
    private final String program = "python";
    private final String options = "-u";

    @Override
    public void createTempNumpyFile(String pythonPath, String uniqueNameNumpy, String pdbPath, String parameter) {
        try {
            ProcessBuilder pb = new ProcessBuilder()
                    .command(program, options, pythonPath, pdbPath, uniqueNameNumpy, parameter);
            Process p = pb.start();
            if (p.waitFor() != EXIT_CODE) {
                LOGGER.warning("There was an error while creating a temporary numpy file");
            }
        } catch (IOException | InterruptedException ex) {
            LOGGER.severe("Error while creating a numpy file, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    public String getChainsPBD(String pythonPath, String pdbID) {
        try {
            ProcessBuilder pb = new ProcessBuilder().command(program, options, pythonPath, pdbID);
            Process p = pb.start();
            BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line;
            StringBuilder buffer = new StringBuilder();
            while ((line = in.readLine()) != null){
                buffer.append(line);;
            }
            if (p.waitFor() != EXIT_CODE) {
                LOGGER.warning("There was an error while retrieving the chains");
            }
            in.close();
            return buffer.toString();
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public String createPcaPlot(String pythonPath, String numpyPath) {
        ProcessBuilder pb = new ProcessBuilder()
                .command(program, options, pythonPath, numpyPath);
        try {
            Process p = pb.start();
            BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            StringBuilder buffer = new StringBuilder();
            String line;
            while ((line = in.readLine()) != null){
                buffer.append(line);
            }
            if (p.waitFor() != 0) {
                LOGGER.warning("There was an error while creating pca plot");
            }
            in.close();
            return buffer.toString();
        } catch (IOException | InterruptedException ex) {
            LOGGER.warning("Error while reading creating pca plot, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    @Override
    public String createScatterPlot(String pythonPath, String numpyPath) {
        ProcessBuilder pb = new ProcessBuilder()
                .command(program, options, pythonPath, numpyPath);
        try {
            Process p = pb.start();
            BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            StringBuilder buffer = new StringBuilder();
            String line;
            while ((line = in.readLine()) != null){
                buffer.append(line);
            }
            if (p.waitFor() != 0) {
                LOGGER.warning("There was an error while creating scatter plot");
            }
            in.close();
            return buffer.toString();
        } catch (IOException | InterruptedException ex) {
            LOGGER.warning("Error while reading creating scatter plot, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    @Override
    public String createPlotlyPcaPlot(String pythonPath, String numpyPath) {
        ProcessBuilder pb = new ProcessBuilder()
                .command(program, options, pythonPath, numpyPath);
        try {
            Process p = pb.start();
            BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            StringBuilder buffer = new StringBuilder();
            String line;
            while ((line = in.readLine()) != null){
                buffer.append(line);
            }
            if (p.waitFor() != 0) {
                LOGGER.warning("There was an error while creating PCA plot");
            }
            in.close();
            return buffer.toString();
        } catch (IOException | InterruptedException ex) {
            LOGGER.warning("Error while reading creating PCA plot, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }


    @Override
    public String PDBAnalyse(String pythonPath, String pdbID, String param, String comparePDB) {
        try {
            ProcessBuilder pb = new ProcessBuilder().command(program, options, pythonPath, pdbID, param, comparePDB);
            Process p = pb.start();
            BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line;
            StringBuilder buffer = new StringBuilder();
            while ((line = in.readLine()) != null){
                buffer.append(line);;
            }
            if (p.waitFor() != 0) {
                LOGGER.warning("There was an error while testing the chains");
            }
            in.close();
            return buffer.toString();
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }
}
