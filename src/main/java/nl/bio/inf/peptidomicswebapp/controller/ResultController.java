package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletRequest;
import jakarta.servlet.http.HttpSession;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.models.PDB;
import nl.bio.inf.peptidomicswebapp.models.Plot;
import nl.bio.inf.peptidomicswebapp.service.PythonService;
import org.springframework.core.io.ClassPathResource;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.RestController;

import java.io.File;
import java.io.IOException;
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

    public ResultController(PythonService pythonService) {
        this.pythonService = pythonService;
    }

    /**
     * This method will create a temp file. And set the file in the session.
     * @param request
     * @param session
     */
    @PostMapping(value = "/create_temp_file")
    public void createTemporaryFile(HttpServletRequest request, HttpSession session) {
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
     */
    @PostMapping(value = "/create_compare_temp")
    public Plot createTemporaryFileCompare(HttpServletRequest request, HttpSession session) {
        try {
            String compareCode = String.valueOf(request.getSession().getAttribute("compareCode"));
            String location = PDB.createTempFile(compareCode);
            session.setAttribute("tempLocationCompare", location);
            File folderScripts = new ClassPathResource("scripts").getFile();
            File fullPath = null;
            // Get the location for the python file
            for (File f: folderScripts.listFiles()) {
                if("PDBAnalyse.py".equals(f.getName())) {
                    fullPath = f;
                }
            }
            // Call the python script
            String bytes = pythonService.PDBAnalyse(
                    fullPath.toString(),
                    request.getSession().getAttribute("tempLocation").toString(),
                    request.getSession().getAttribute("parameter").toString(),
                    location
            );
            return new Plot(bytes);
        } catch (IOException ex) {
            LOGGER.warning("Error while performing the script on the data, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }

    }

    /**
     * This method will call a python script that will retrieve the chains of the pdb file and return to the site.
     * @param request
     * @return
     */
    @PostMapping(value = "/get_chains")
    public @ResponseBody Plot getChains(HttpServletRequest request){
        try {
            File folderScripts  = new ClassPathResource("scripts").getFile();
            File fullPath = null;
            // Get the location for the python file
            for (File f: folderScripts.listFiles()) {
                if("retrieve_chains_pdb.py".equals(f.getName())) {
                    fullPath = f;
                }
            }
            // Call the python script
            String chain = pythonService.getChainsPBD(
                    String.valueOf(fullPath),
                    request.getSession().getAttribute("tempLocation").toString());
            return new Plot(chain);
        } catch (IOException ex) {
            LOGGER.warning("Error while retrieving the chains from the data, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }
}