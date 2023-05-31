package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletRequest;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.exceptions.TooLargeNumberException;
import org.apache.tomcat.util.http.fileupload.impl.SizeLimitExceededException;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.ControllerAdvice;
import org.springframework.web.bind.annotation.ExceptionHandler;

import java.util.logging.Logger;

/**
 * @author Wouter Zeevat
 */

@ControllerAdvice
public class SizeLimitExceptionController {

    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());

    /**
     * This method redirects the user to the pdb error page when a sizeLimitExceedException occurs
     *
     * @param model
     * @return
     */
    @ExceptionHandler(SizeLimitExceededException.class)
    public String handleMaxFile(Model model, HttpServletRequest request) {
        LOGGER.severe(request.getSession().getId()+ " session id tried to upload a file that is too large");
        model.addAttribute("code", "500");
        model.addAttribute("message", "This file is too large for the program to handle!");
        return "pdb_error";
    }

    @ExceptionHandler(TooLargeNumberException.class)
    public String handleMaxOligo(Model model) {
        model.addAttribute("code", "500");
        model.addAttribute("message", "This oligo number is too large!");
        return "pdb_error";
    }

}
