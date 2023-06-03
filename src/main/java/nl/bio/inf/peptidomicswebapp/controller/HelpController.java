package nl.bio.inf.peptidomicswebapp.controller;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;

/**
 * Controller that handles the tutorial web page
 *
 * @author Wouter Zeevat
 */
@Controller
public class HelpController {

    @Value("${spring.servlet.multipart.max-request-size}")
    private String maxMB;

    @Value("${max.oligo.length}")
    private int maxOligo;

    @RequestMapping(value = {"/help"})
    public String landingPage(Model model){
        model.addAttribute("oligoLength", maxOligo);
        model.addAttribute("maxMB", maxMB);
        return "help";
    }
}
