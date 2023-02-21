package nl.bio.inf.peptidomicswebapp.controller;

import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.GetMapping;


@Controller
public class UploadController {
    @GetMapping(value ="/upload")
    public String landingPage(){
        return "upload";
    }
}
