package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.RequestDispatcher;
import jakarta.servlet.http.HttpServletRequest;
import nl.bio.inf.peptidomicswebapp.service.ErrorService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.boot.web.servlet.error.ErrorController;

@Controller
public class CustomErrorController implements ErrorController {
    @Autowired
    private ErrorService errorService;

    @RequestMapping(value = "/error")
    public String renderErrorPage(Model model, final HttpServletRequest request) {
        final int errorCode = getHttpStatusCode(request);
        final String errorMessage = errorService.generateErrorMessage(errorCode);
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
