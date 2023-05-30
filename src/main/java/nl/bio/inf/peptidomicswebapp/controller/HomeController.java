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
public class HomeController {

    private static final String[] WORKFLOWLIST = new String[]{
            "Retrieve the coordinates and other information like peptide from the PDB file",
            "Use the XSSP api to get the secondary structures",
            "Split the coordinates on breaks/chains",
            "Get the internal distance between the coordinates",
            "Get the cov matrix of the distances",
            "Perform PCA on the cov matrix and retrieve the eigenvectors",
            "Perform matrix multiplication on the distance and eigenvectors"
    };

    @RequestMapping(value = {"", "/", "/home"})
    public String landingPage(Model model){
        model.addAttribute("workflowList", WORKFLOWLIST);
        return "index";
    }
}