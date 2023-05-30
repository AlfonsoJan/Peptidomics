package nl.bio.inf.peptidomicswebapp.controller;

import org.apache.tomcat.util.http.fileupload.impl.SizeLimitExceededException;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.ControllerAdvice;
import org.springframework.web.bind.annotation.ExceptionHandler;

/**
 * @author Wouter Zeevat
 */

@ControllerAdvice
public class SizeLimitExceptionController {

    /**
     * This method redirects the user to the pdb error page when a sizeLimitExceedException occurs
     *
     * @param model
     * @return
     */
    @ExceptionHandler(SizeLimitExceededException.class)
    public String handleMaxFile(Model model) {
        model.addAttribute("code", "500");
        model.addAttribute("message", "This file is too large for the program to handle!");
        return "pdb_error";
    }

}
