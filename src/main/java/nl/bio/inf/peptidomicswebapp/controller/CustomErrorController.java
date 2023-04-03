package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.RequestDispatcher;
import jakarta.servlet.http.HttpServletRequest;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.config.HtmlLogFormatter;
import nl.bio.inf.peptidomicswebapp.service.ErrorService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.boot.web.servlet.error.ErrorController;

import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;

@Controller
public class CustomErrorController implements ErrorController {

    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());
    @Autowired
    private ErrorService errorService;

    @RequestMapping(value = "/error")
    public String renderErrorPage(Model model, final HttpServletRequest request) throws IOException {
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

    private int getHttpStatusCode(final HttpServletRequest request) {
        if (request.getAttribute(RequestDispatcher.ERROR_STATUS_CODE) == null) {
            return -1;
        }
        return (int) request.getAttribute(RequestDispatcher.ERROR_STATUS_CODE);
    }
}