package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletResponse;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import org.springframework.web.servlet.ModelAndView;

import java.io.IOException;
import java.util.Arrays;


@Controller
public class UploadController {
    @GetMapping(value ="/upload")
    public String landingPage(){
        return "upload";
    }

    @PostMapping(value = "/result_from_code")
    public String resultFromCode(@RequestParam("pdb_code") String text) {
        System.out.println(text);
        return "upload";
    }
    @PostMapping(value = "/result_from_files")
    public String resultFromFiles(@RequestParam("pdb_file") MultipartFile file) {

//        Arrays.asList(files).stream().forEach(file -> System.out.println(file.getOriginalFilename()));
//        try {
//            BufferedReader reader = new BufferedReader(new InputStreamReader(file.getInputStream()));
//            String strCurrentLine;
//            while ((strCurrentLine = reader.readLine()) != null) {
//                System.out.println(strCurrentLine);
//            }
//        } catch (IOException e) {
//            throw new RuntimeException(e);
//        }
        return "upload";
    }
}
