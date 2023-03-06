package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletResponse;
import jakarta.servlet.http.HttpSession;
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
    public String resultFromCode(@RequestParam("pdb_code") String text) {
        System.out.println(text);
        return "redirect:/result";
    }
    @PostMapping(value = "/result_from_files")
    public String resultFromFiles(HttpSession session,
                                  @RequestParam("pdb_file") MultipartFile file,
                                  Model model) throws IOException {
        model.asMap().put("uploadedFile", new UploadedFile(file.getOriginalFilename(), file.getBytes()));
        return "redirect:/result";
    }

    @RequestMapping(value = "/result")
    public String resultFromFiles(Model model) {
        try {
            UploadedFile uploadedFile = (UploadedFile) model.getAttribute("uploadedFile");
            System.out.println(uploadedFile.originalFilename());
            return "index";
        }
        catch (ClassCastException ex) {
            return "redirect:/";
        }
    }
}
