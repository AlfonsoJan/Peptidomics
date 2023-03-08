package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletResponse;
import jakarta.servlet.http.HttpSession;
import nl.bio.inf.peptidomicswebapp.models.PDB;
import nl.bio.inf.peptidomicswebapp.models.UploadedFile;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import org.springframework.web.servlet.ModelAndView;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


@Controller
@SessionAttributes(value = {"uploadedFile"})
public class UploadController {
    @GetMapping(value ="/upload")
    public String landingPage(){
        return "upload";
    }

    @PostMapping(value = "/result_from_code")
    public String resultFromCode(@RequestParam("pdb_code") String pdbCode, Model model) throws IOException {
        model.asMap().put("uploadedFile", "");
        PDB pdb = new PDB(pdbCode);
        String filename = pdb.getStructureId() + ".pdb";
        model.asMap().put("uploadedFile", new UploadedFile(filename, pdb.getBytes()));
        return "redirect:/result";
    }

    @PostMapping(value = "/result_from_files")
    public String resultFromFiles(@RequestParam("pdb_file") MultipartFile file,
                                  Model model) throws IOException {
        model.asMap().put("uploadedFile", new UploadedFile(file.getOriginalFilename(), file.getBytes()));
        return "redirect:/result";
    }

    @RequestMapping(value = "/result")
    public String resultPage(Model model) {
        try {
            UploadedFile uploadedFile = (UploadedFile) model.getAttribute("uploadedFile");
            if (uploadedFile != null) {
                System.out.println(uploadedFile.getBytes());
                return "results";
            }
            return "redirect:/";
        }
        catch (ClassCastException ex) {
            return "redirect:/";
        }
    }
}
