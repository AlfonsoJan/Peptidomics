package nl.bio.inf.peptidomicswebapp.controller;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;

/**
 *  This class handles request to the home page.
 * @author Jan Alfonso Busker
 */
@Controller
public class HomeController {

    @Value("${home.page.project.info.text}")
    private String infoText;

    @RequestMapping(value = {"", "/", "/home"})
    public String landingPage(Model model){
        model.addAttribute("infoText", infoText);
        return "index";
    }
}