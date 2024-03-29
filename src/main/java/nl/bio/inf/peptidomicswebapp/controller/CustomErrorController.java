package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.RequestDispatcher;
import jakarta.servlet.http.HttpServletRequest;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.service.ErrorService;
import org.springframework.boot.web.servlet.error.ErrorController;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;

import java.util.logging.Logger;

/**
 *  This class handles the errors and what message to show on the page.
 * @author Wouter Zeevat
 */

@Controller
public class CustomErrorController implements ErrorController {

    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());
    private final ErrorService errorService;

    public CustomErrorController(ErrorService errorService) {
        this.errorService = errorService;
    }

    /**
     * Renders the error page with the errorCode and message according to the errorService
     * @param model model
     * @param request request
     * @return String name for the html page
     */
    @RequestMapping(value = "/error")
    public String renderErrorPage(Model model, final HttpServletRequest request) {
        final int errorCode = getHttpStatusCode(request);
        final String errorMessage = errorService.generateErrorMessage(errorCode);
        String originalUrl = (String) request.getAttribute(RequestDispatcher.FORWARD_REQUEST_URI);
        if (404 == errorCode) {
            LOGGER.severe(String.format("Error with code: 404, url: %s", originalUrl));
        } else {
            LOGGER.severe(String.format("Error with code: %s", errorCode));
        }
        model.addAttribute("errorCode", errorCode);
        model.addAttribute("errorMsg", errorMessage);
        return "error";
    }

    /**
     * Renders the custom PDB error page with the errorCode and message according to the errorService
     * @param model model
     * @param request request
     * @return String for the html page
     */
    @RequestMapping(value = "/pdb_error")
    public String renderPDBErrorPage(Model model, final HttpServletRequest request) {
        final String error = request.getParameter("code");
        String message = request.getParameter("message");

        model.addAttribute("code", error);
        model.addAttribute("message", message);
        return "pdb_error";
    }

    /**
     * Returns the error status code
     * @param request request
     * @return integer
     */
    private int getHttpStatusCode(final HttpServletRequest request) {
        if (request.getAttribute(RequestDispatcher.ERROR_STATUS_CODE) == null) {
            return -1;
        }
        return (int) request.getAttribute(RequestDispatcher.ERROR_STATUS_CODE);
    }
}