package nl.bio.inf.peptidomicswebapp.service;

import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Service;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.logging.Logger;

@Service
public class PythonService implements PythonConstructor{

    @Value("${python.path-name}")
    private String program;

    @Value("${python.path-optons}")
    private String options;
    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());
    private static final int EXIT_CODE = 0;

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