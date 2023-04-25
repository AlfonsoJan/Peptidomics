package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletRequest;
import jakarta.servlet.http.HttpSession;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.models.PDB;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;

import java.io.IOException;
import java.util.logging.Logger;

/**
 *  This class handles the file/code uploads. And set these in a session
 * @author Jan Alfonso
 * @author Wouter Zeevat
 */

@Controller
public class UploadController {

    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());
    @GetMapping(value ="/upload")
    public String landingPage(){
        return "upload";
    }

    /**
     * This method creates a session and sets the codes in the session and redirect to the result page.
     * @param pdbCode
     * @param paramCode
     * @param compareCode
     * @param session
     */
    @PostMapping(value = "/result_from_code")
    public String resultFromCode(@RequestParam("pdb_code") String pdbCode,
                                 String paramCode,
                                 String compareCode,
                                 HttpSession session) {
        try {
            PDB testPDB = new PDB(pdbCode);
            session.setAttribute("parameter", paramCode);
            session.setAttribute("PDBFiles", testPDB);
            session.setAttribute("compareCode", compareCode);
            return "redirect:/result";
        } catch (IOException ex) {
            LOGGER.warning("Error while reading creating PDB class with pdb code, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    /**
     * This method creates a session and sets the file and codes in the session and redirect to the result page.
     * @param file
     * @param paramFile
     * @param compareFile
     * @param session
     */
    @PostMapping(value = "/result_from_files")
    public String resultFromFiles(@RequestParam("pdb_file") MultipartFile file,
                                  String paramFile,
                                  String compareFile,
                                  HttpSession session) {
        try {
            PDB pdb = new PDB(file.getBytes(), file.getOriginalFilename());
            session.setAttribute("parameter", paramFile);
            session.setAttribute("PDBFiles", pdb);
            session.setAttribute("compareCode", compareFile);
            return "redirect:/result";
        } catch (IOException ex) {
            LOGGER.warning("Error while reading PDB file, message=" + ex.getMessage());
            throw new RuntimeException(ex);

        }
    }

    /**
     * This method will return the result page if its a correct pdb file/code.
     * And if there is nothing in the session then go to the upload page
     * @param model
     * @param request
     */
    @RequestMapping(value = "/result")
    public String resultPage(Model model, HttpServletRequest request){
        try {
            // If the session if null, then redirect to the upload page
            PDB pdb = (PDB) request.getSession().getAttribute("PDBFiles");
            if (pdb == null || pdb.getStructureId() == null) {
                LOGGER.warning(String.format("PDB structure of %s is null", request.getSession().getId()));
                return "redirect:/upload";
            }
            model.addAttribute("fileName", "<strong>Results of: </strong>" + pdb.getStructureId());
            return "results";
        } catch (ClassCastException ex) {
            LOGGER.warning("Error while class casting to PDB, message=" + ex.getMessage());
            return "redirect:/upload";
        }
    }
}