package nl.bio.inf.peptidomicswebapp.controller;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.PropertySource;
import org.springframework.core.env.Environment;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;

/**
 *  This class handles request to the home page.
 * @author Jan Alfonso Busker
 */
@Controller
@PropertySource("classpath:index.properties")
public class HomeController {

    @Autowired
    private Environment env;


    @RequestMapping(value = {"", "/", "/home"})
    public String landingPage(Model model){
        model.addAttribute("infoTxt", env.getProperty("info.txt"));
        return "index";
    }
}