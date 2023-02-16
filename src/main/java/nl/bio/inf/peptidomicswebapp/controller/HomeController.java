package nl.bio.inf.peptidomicswebapp.controller;

import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;

@Controller
public class HomeController {
    @RequestMapping(value = {"", "/", "/home"})
    public String landingPage(Model model){
        model.addAttribute("brand", "Peptidomics");
        model.addAttribute("title", "Peptidomics");
        return "index";
    }
}
