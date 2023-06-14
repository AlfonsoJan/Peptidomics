package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletRequest;
import jakarta.servlet.http.HttpSession;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.exceptions.InvalidPDBCodeException;
import nl.bio.inf.peptidomicswebapp.exceptions.TooLargeNumberException;
import nl.bio.inf.peptidomicswebapp.models.PDB;
import org.apache.tomcat.util.http.fileupload.impl.SizeLimitExceededException;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.multipart.MultipartFile;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.logging.Logger;

/**
 *  This class handles the file/code uploads. And set these in a session
 * @author Jan Alfonso
 * @author Wouter Zeevat
 */

@Controller
public class UploadController {

    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());

    @Value("${spring.servlet.multipart.max-request-size}")
    private String maxMB;

    @Value("${max.oligo.length}")
    private int maxOligo;


    @GetMapping(value ="/upload")
    public String landingPage(Model model){
        model.addAttribute("maxOligo", maxOligo);
        return "upload";
    }

    /**
     * This method creates a session and sets the codes in the session and redirect to the result page.
     *
     * @param pdbCode pdb code
     * @param oligoParam oligo length
     * @param session session
     * @throws RuntimeException when the PDB is invalid
     */
    @PostMapping(value = "/result_from_code")
    public String resultFromCode(@RequestParam("pdb_code") String pdbCode,
                                 String oligoParam,
                                 HttpSession session) {

        try {
            PDB pdb = new PDB(pdbCode);
            deletePreviousSession(session);

            if (!pdb.isValid()) {
                return ("redirect:/pdb_error?code=" + pdb.getFileName());
            }

            if (Integer.parseInt(oligoParam) > maxOligo) throw new TooLargeNumberException();

            // If pdb file size is higher than the max amount
            if ((pdb.getBytes().length) / (1024 * 1024) > Integer.parseInt(maxMB.toLowerCase().replace("mb", ""))) {
                LOGGER.severe("File too large!");
                throw new SizeLimitExceededException("pdb code (" + pdb.getStructureId() + ") too large!", (pdb.getBytes().length) / (1024 * 1024), 10);
            }

            session.setAttribute("pepSize", oligoParam);
            session.setAttribute("PDBFiles", pdb);
            return "redirect:/result";
        } catch (IOException | InvalidPDBCodeException ex) {
            LOGGER.warning("Error while reading creating PDB class with pdb code, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        } catch (TooLargeNumberException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * This method creates a session and sets the file and codes in the session and redirect to the result page.
     * @param file uploaded file
     * @param oligoParam oligo length
     * @param session session
     * @throws RuntimeException when PDB file can not be read correctly
     */
    @PostMapping(value = "/result_from_files")
    public String resultFromFiles(@RequestParam("pdb-file") MultipartFile file,
                                  String oligoParam,
                                  HttpSession session) {
        try {
            // Creates PDB instance and redirects to page
            PDB pdb = new PDB(file.getBytes(), file.getOriginalFilename());
            deletePreviousSession(session);

            if (!pdb.isValid()) {
                return ("redirect:/pdb_error?code=1f6e");
            }

            if (Integer.parseInt(oligoParam) > maxOligo) throw new TooLargeNumberException();

            session.setAttribute("pepSize", oligoParam);
            session.setAttribute("PDBFiles", pdb);
            return "redirect:/result";
        } catch (IOException ex) {
            LOGGER.warning("Error while reading PDB file, message=" + ex.getMessage());
            throw new RuntimeException(ex);

        } catch (TooLargeNumberException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * This method will return the result page if it's a correct pdb file/code.
     * And if there is nothing in the session then go to the upload page
     * @param model model
     * @param request request
     * @throws ClassCastException when PDB can't turn into a PDB class instance
     */
    @RequestMapping(value = "/result")
    public String resultPage(Model model, HttpServletRequest request){
        try {
            // If the session is null, then redirect to the upload page
            PDB pdb = (PDB) request.getSession().getAttribute("PDBFiles");
            if (pdb == null || pdb.getStructureId() == null) {
                LOGGER.warning(String.format("PDB structure of %s is null", request.getSession().getId()));
                return "redirect:/upload";
            }
            model.addAttribute("fileName", pdb.getStructureId());
            model.addAttribute("pepSize", request.getSession().getAttribute("pepSize"));
            return "results";
        } catch (ClassCastException ex) {
            LOGGER.warning("Error while class casting to PDB, message=" + ex.getMessage());
            return "redirect:/upload";
        }
    }

    /**
     * This method checks if there is an existing session and deletes it if there is one.
     *
     * @param session session
     * @throws RuntimeException when the PDB is invalid
     */
    public void deletePreviousSession(HttpSession session) {
        if (session.getAttribute("tempLocation") != null) {
            String tempLocation = String.valueOf(session.getAttribute("tempLocation"));
            String jsonFileLocation = String.valueOf(session.getAttribute("jsonFile"));
            try {
                Files.delete(Path.of(tempLocation));
                Files.delete(Path.of(jsonFileLocation));
                LOGGER.info("Deleted session files of: " + session.getId());
            } catch (IOException e) {
                LOGGER.warning("Could not delete the files of session: " + session.getId());
                throw new RuntimeException(e);
            }
            session.removeAttribute("chains");
            session.removeAttribute("tempLocation");
            session.removeAttribute("analysis");
        }

    }
}