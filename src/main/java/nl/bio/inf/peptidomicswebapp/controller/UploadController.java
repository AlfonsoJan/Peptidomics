package nl.bio.inf.peptidomicswebapp.controller;

import nl.bio.inf.peptidomicswebapp.models.PDB;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;

import java.io.IOException;


@Controller
@SessionAttributes(value = {"PDBFiles"})
public class UploadController {
    @GetMapping(value ="/upload")
    public String landingPage(){
        return "upload";
    }

    @PostMapping(value = "/result_from_code")
    public String resultFromCode(@RequestParam("pdb_code") String pdbCode, Model model) throws IOException {
        model.asMap().put("PDBFiles", new PDB(pdbCode));
        return "redirect:/result";
    }

    @PostMapping(value = "/result_from_files")
    public String resultFromFiles(@RequestParam("pdb_file") MultipartFile file,
                                  Model model) throws IOException {
        model.asMap().put("PDBFiles", new PDB(file.getOriginalFilename(), file.getBytes(), PDB.getStructureFromInputstream(file.getInputStream())));
        return "redirect:/result";
    }

    @RequestMapping(value = "/result")
    public String resultPage(Model model) {
        try {
            PDB pdbFile = (PDB) model.getAttribute("PDBFiles");
            if (pdbFile != null) {
                String fileStructure = pdbFile.getStructureId();
                if (pdbFile.getStructureId() ==  null) {
                    fileStructure = "''";
                }
                model.addAttribute("fileName", "<strong class=\"is-size-2\">Results of: </strong>" + fileStructure);
                return "results";
            }
            return "redirect:/";
        }
        catch (ClassCastException ex) {
            return "redirect:/";
        }
    }
}
