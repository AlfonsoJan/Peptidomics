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


@Controller
public class UploadController {

    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());
    @GetMapping(value ="/upload")
    public String landingPage(){
        return "upload";
    }

    @PostMapping(value = "/result_from_code")
    public String resultFromCode(@RequestParam("pdb_code") String pdbCode,
                                 String param_code,
                                 String compare_code,
                                 HttpSession session) {
        try {
            PDB pdb = new PDB(pdbCode);
            session.setAttribute("parameter", param_code);
            session.setAttribute("PDBFiles", pdb);
            session.setAttribute("compareCode", compare_code);
            return "redirect:/result";
        } catch (IOException ex) {
            LOGGER.warning("Error while reading creating PDB class with pdb code, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    @PostMapping(value = "/result_from_files")
    public String resultFromFiles(@RequestParam("pdb_file") MultipartFile file,
                                  String param_file,
                                  String compare_file,
                                  HttpSession session) {
        session.setAttribute("compareCode", compare_file);
        try {
            PDB pdb = new PDB(file.getOriginalFilename(), file.getBytes(), PDB.getStructureFromInputstream(file.getInputStream()));
            session.setAttribute("parameter", param_file);
            session.setAttribute("PDBFiles", pdb);
            return "redirect:/result";
        } catch (IOException ex) {
            LOGGER.warning("Error while reading PDB file, message=" + ex.getMessage());
            throw new RuntimeException(ex);

        }
    }

    @RequestMapping(value = "/result")
    public String resultPage(Model model, HttpServletRequest request){
        try {
            PDB pdb = (PDB) request.getSession().getAttribute("PDBFiles");
            if (pdb == null || pdb.getStructureId() == null) {
                LOGGER.warning(String.format("PDB structure of %s is null", request.getSession().getId()));
                return "redirect:/";
            }
            model.addAttribute("fileName", "<strong>Results of: </strong>" + pdb.getStructureId());
            return "results";
        } catch (ClassCastException ex) {
            LOGGER.warning("Error while class casting to PDB, message=" + ex.getMessage());
            return "redirect:/";
        }
    }
}