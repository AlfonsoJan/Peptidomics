package nl.bio.inf.peptidomicswebapp.controller;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.PropertySource;
import org.springframework.core.env.Environment;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;

import java.util.Arrays;

/**
 *  This class handles request to the home page.
 * @author Jan Alfonso Busker
 */
@PropertySource("classpath:messages_en.properties")
@Controller
public class HomeController {

    @Value("${workflow.txt}")
    private String workflowText;

    @Value("${pdb.files}")
    private String pdbFiles;

    @RequestMapping(value = {"", "/", "/home"})
    public String landingPage(Model model){
        model.addAttribute("workflowList", workflowText.split(";"));
        model.addAttribute("eigenVectors", pdbFiles.split(";"));
        return "index";
    }
}