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

    private static final String[] eigenVectors = new String[]{"6cw5", "7juk", "3qo2", "4ylk", "2wef", "3v1f", "7xoc", "6njn", "3w00", "6fn1",
            "6eth", "3u9a", "8aop", "2ppv", "4qhs", "3jr7", "4keq", "2x93", "5umf", "1nv1", "6ivy", "4pjm", "3v7b", "7sb5", "2fur", "5xur",
            "6utd", "7xix", "5z4o", "7k6w", "6vhh", "5pct", "5x1m", "4cug", "2qt6", "1ck6", "1igf",
            "5pt0", "7e7c", "3o0f", "1bcr", "6mri", "7koz", "6mr5",
            "6uki", "6bpp", "1zxi", "4dm0", "5w7e", "6rgj"};

    @RequestMapping(value = {"", "/", "/home"})
    public String landingPage(Model model){
        model.addAttribute("workflowList", WORKFLOWLIST);
        model.addAttribute("eigenVectors", eigenVectors);
        return "index";
    }
}