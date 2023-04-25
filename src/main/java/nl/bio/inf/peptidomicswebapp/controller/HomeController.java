package nl.bio.inf.peptidomicswebapp.controller;

import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestMapping;

/**
 *  This class handles request to the home page.
 * @author Jan Alfonso
 */
@Controller
public class HomeController {
    @RequestMapping(value = {"", "/", "/home"})
    public String landingPage(){
        return "index";
    }
}