package nl.bio.inf.peptidomicswebapp.service;

import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.exceptions.InvalidPDBCodeException;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Service;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.logging.Logger;

/**
 * This class will call the python script and return the scripts.
 * @author Wouter Zeevat
 * @author Jan Alfonso Busker
 */
@Service
public class PythonService implements PythonConstructor{

    @Value("${python.path-name}")
    private String program;

    @Value("${python.path-optons}")
    private String options;
    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());
    private static final int EXIT_CODE = 0;

    /**
     * Calls the python script to get the chains of the pdb file
     * @param pythonPath
     * @param pdbPath
     * @return
     */
    @Override
    public String getChainsPBD(String pythonPath, String pdbPath) {
        try {
            if (program == null) program = "python3";
            if (options == null) options = "-u";

            if (pythonPath == null | pdbPath == null) throw new NullPointerException();

            ProcessBuilder pb = new ProcessBuilder().command(program, options, pythonPath, pdbPath);
            Process p = pb.start();
            BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line;
            StringBuilder buffer = new StringBuilder();
            while ((line = in.readLine()) != null){
                buffer.append(line);;
            }
            if (p.waitFor() != EXIT_CODE) {
                LOGGER.warning("There was an error while retrieving the chains");
                throw new FileNotFoundException();
            }
            in.close();
            return buffer.toString();
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Calls the python script to perform PCA analysis and retrieve the result
     * @param pythonPath
     * @param pdbPath
     * @return
     */
    @Override
    public String PDBAnalyse(String pythonPath, String pdbPath, String pepSize, String pdbCode) throws InvalidPDBCodeException {
        try {
            if (program == null) program = "python3";
            if (options == null) options = "-u";

            if (pythonPath == null | pdbPath == null | pepSize == null | pdbCode == null) throw new NullPointerException();

            int size = Integer.parseInt(pepSize);
            if (size < 1 | size > 30) throw new NumberFormatException();
            if (pdbCode.length() != 4) throw new InvalidPDBCodeException();

            ProcessBuilder pb = new ProcessBuilder().command(program, options, pythonPath, pdbPath, pepSize, pdbCode);
            Process p = pb.start();
            BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line;
            StringBuilder buffer = new StringBuilder();
            while ((line = in.readLine()) != null){
                buffer.append(line);;
            }
            if (p.waitFor() != 0) {
                LOGGER.warning("There was an error while performing the analysis!");
                throw new FileNotFoundException();
            }
            in.close();
            return buffer.toString();
        } catch (IOException | InterruptedException | InvalidPDBCodeException e){
            throw new InvalidPDBCodeException();
        }
    }
}