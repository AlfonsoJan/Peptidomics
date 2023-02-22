package nl.bio.inf.peptidomicswebapp.controller;

import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestMapping;

@Controller
public class UploadController {
    @RequestMapping(value = {"/results"})
    public String landingPage(){
        return "results";
    }
}
