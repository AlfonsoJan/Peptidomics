package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletRequest;
import jakarta.servlet.http.HttpSession;
import nl.bio.inf.peptidomicswebapp.models.PDB;
import org.springframework.core.io.ClassPathResource;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;


@Controller
public class UploadController {
    @GetMapping(value ="/upload")
    public String landingPage(){
        return "upload";
    }

    @PostMapping(value = "/result_from_code")
    public String resultFromCode(@RequestParam("pdb_code") String pdbCode,
                                 String param_code,
                                 HttpSession session) throws IOException {
        session.setAttribute("parameter", param_code);
        PDB pdb = new PDB(pdbCode);
        session.setAttribute("PDBFiles", pdb);
        return "redirect:/result";
    }

    @PostMapping(value = "/result_from_files")
    public String resultFromFiles(@RequestParam("pdb_file") MultipartFile file,
                                  String param_file,
                                  HttpSession session) throws IOException {
        PDB pdb = new PDB(file.getOriginalFilename(), file.getBytes(), PDB.getStructureFromInputstream(file.getInputStream()));
        session.setAttribute("parameter", param_file);
        session.setAttribute("PDBFiles", pdb);
        return "redirect:/result";
    }

    @RequestMapping(value = "/result")
    public String resultPage(Model model, HttpServletRequest request) throws IOException {
        File folderScripts = new ClassPathResource("scripts").getFile();
        try {
            PDB pdb = (PDB) request.getSession().getAttribute("PDBFiles");
            if (pdb == null || pdb.getStructureId() == null) {
                return "redirect:/";
            }
            model.addAttribute("fileName", "<strong>Results of: </strong>" + pdb.getStructureId());
            return "results";
        } catch (ClassCastException ex) {
            return "redirect:/";
        }
    }
}