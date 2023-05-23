package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletRequest;
import jakarta.servlet.http.HttpSession;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.exceptions.InvalidPDBCodeException;
import nl.bio.inf.peptidomicswebapp.models.PDB;
import nl.bio.inf.peptidomicswebapp.models.Plot;
import nl.bio.inf.peptidomicswebapp.service.PythonService;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.springframework.security.web.csrf.CsrfToken;
import org.springframework.web.bind.annotation.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Logger;

/**
 *  This class handles the result methods.
 * @author Jan Alfonso Busker
 * @author Wouter Zeevat
 */

@RestController
public class ResultController {
    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());

    private final PythonService pythonService;
    @Value("${python.executable.folder}")
    private String pythonFolder;

    public ResultController(PythonService pythonService) {
        this.pythonService = pythonService;
    }


    /**
     * This will return the csrf token for the stateless fetch calls
     * @param request
     * @return csrf token
     */
    @RequestMapping(value="/csrf-token", method= RequestMethod.GET)
    public @ResponseBody String getCsrfToken(HttpServletRequest request) {
        CsrfToken token = (CsrfToken)request.getAttribute(CsrfToken.class.getName());
        return token.getToken();
    }

    /**
     * This method will create a temp file. And set the file in the session.
     * @param request
     * @param session
     * @throws Exception when file cant be created
     */
    @PostMapping(value = "/create_temp_file")
    public void createTemporaryFile(HttpServletRequest request, HttpSession session) {
        if (session.getAttribute("tempLocation") != null) {
            return;
        }
        try {
            PDB pdb = (PDB) request.getSession().getAttribute("PDBFiles");
            String tempLocation = pdb.createTempFile();
            session.setAttribute("tempLocation", tempLocation);
        } catch (Exception ex) {
            LOGGER.severe("Error while creating a temp file, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    /**
     * This method will create a temp file of the compare pdb code and
     * call the script that will run the analysis and return to the site.
     * @param request
     * @param session
     * @return
     * @throws IOException when the script can't be run correctly
     */
    @PostMapping(value = "/perform_pca_analysis")
    public @ResponseBody Plot performPCAAnalysis(HttpServletRequest request, HttpSession session) {
        if (request.getSession().getAttribute("analysis") != null) {
            return (Plot) session.getAttribute("analysis");
        }
        try {
            PDB pdb = (PDB) request.getSession().getAttribute("PDBFiles");
            Path filePath = Paths.get(pythonFolder, "scripts", "pdb_analysis.py");
            // Call the python script
            String bytes = pythonService.PDBAnalyse(
                    String.valueOf(filePath),
                    request.getSession().getAttribute("tempLocation").toString(),
                    request.getSession().getAttribute("pepSize").toString(),
                    pdb.getStructureId()
            );
            Plot plot = new Plot(bytes);
            request.getSession().setAttribute("analysis", plot);
            return plot;
        } catch (InvalidPDBCodeException ex) {
            LOGGER.warning("Error while performing the script on the data, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    /**
     * This method will call a python script that will retrieve the chains of the pdb file and return to the site.
     * @param request
     * @return
     * @throws RuntimeException when the file can't be read
     */
    @PostMapping(value = "/get_chains")
    public @ResponseBody Plot getChains(HttpServletRequest request){
        if (request.getSession().getAttribute("chains") != null) {
            return (Plot) request.getSession().getAttribute("chains");
        }
        Path filePath = Paths.get(pythonFolder, "scripts", "retrieve_chains_pdb.py");
        // Call the python script
        String chain = pythonService.getChainsPBD(
                String.valueOf(filePath),
                request.getSession().getAttribute("tempLocation").toString());
        Plot plot = new Plot(chain);
        request.getSession().setAttribute("chains", plot);
        return plot;
    }
}