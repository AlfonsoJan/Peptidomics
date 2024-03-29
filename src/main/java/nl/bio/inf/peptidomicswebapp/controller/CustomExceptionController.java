package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletRequest;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.exceptions.EigenVectorsNotFoundException;
import nl.bio.inf.peptidomicswebapp.exceptions.TooLargeNumberException;
import org.apache.tomcat.util.http.fileupload.impl.SizeLimitExceededException;
import org.springframework.http.HttpStatus;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.ControllerAdvice;
import org.springframework.web.bind.annotation.ExceptionHandler;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.ResponseStatus;
import org.springframework.web.method.annotation.MethodArgumentTypeMismatchException;

import java.util.logging.Logger;

/**
 * This class gives advice for the custom exceptions
 * @author Wouter Zeevat
 */

@ControllerAdvice
public class CustomExceptionController {

    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());

    /**
     * This method redirects the user to the pdb error page when a sizeLimitExceedException occurs
     *
     * @param model model
     * @return String for the html page
     */
    @ExceptionHandler(SizeLimitExceededException.class)
    public String handleMaxFile(Model model, HttpServletRequest request) {
        LOGGER.severe(request.getSession().getId()+ " session id tried to upload a file that is too large");
        model.addAttribute("code", "500");
        model.addAttribute("message", "This file is too large for the program to handle!");
        return "pdb_error";
    }

    /**
     * This class handles too large files
     * @param model model
     * @return String for the html page
     */
    @ExceptionHandler(TooLargeNumberException.class)
    public String handleMaxOligo(Model model) {
        model.addAttribute("code", "500");
        model.addAttribute("message", "This oligo number is too large!");
        return "pdb_error";
    }

    /**
     * Handles eigenvectors not found errors for the api
     * @param ex exception
     * @return ErrorResponse
     */
    @ResponseBody
    @ExceptionHandler(EigenVectorsNotFoundException.class)
    @ResponseStatus(HttpStatus.NOT_FOUND)
    ErrorResponse employeeNotFoundHandler(EigenVectorsNotFoundException ex) {
        return new ErrorResponse(HttpStatus.NOT_FOUND.value(), ex.getMessage());
    }

    /**
     * Handles when the user give a string and not a integer for the api
     * @return ErrorResponse
     */
    @ResponseBody
    @ExceptionHandler(MethodArgumentTypeMismatchException.class)
    @ResponseStatus(HttpStatus.BAD_REQUEST)
    ErrorResponse invalidTypeApiHandles() {
        return new ErrorResponse(HttpStatus.BAD_REQUEST.value(), "The length needs to be an integer!");
    }

    /**
     * Class for the error response
     * @param status
     * @param text
     */
    public record ErrorResponse(int status, String text) {}
}
