package nl.bio.inf.peptidomicswebapp.service;

import org.springframework.context.annotation.PropertySource;
import org.springframework.core.env.Environment;
import org.springframework.stereotype.Service;

/**
 * This class will handle all the errors thrown.
 * @author Jan Alfonso Busker
 */
@Service
@PropertySource("classpath:httpErrorCodes.properties")
public class ErrorService {

    private final Environment env;

    public ErrorService(Environment env) {
        this.env = env;
    }

    /**
     * Generates the error text
     * @param errorCode the error code
     * @return String of the error message
     */
    public String generateErrorMessage(final int errorCode) {
        String message = "There was a code " + errorCode + " error: ";
        switch (errorCode) {
            case 302 -> message += env.getProperty("302");
            case 304 -> message += env.getProperty("304");
            case 400 -> message += env.getProperty("400");
            case 401 -> message += env.getProperty("401");
            case 404 -> message += env.getProperty("404");
            case 405 -> message += env.getProperty("405");
            case 410 -> message += env.getProperty("410");
            case 500 -> message += env.getProperty("500");
            case 502 -> message += env.getProperty("502");
            default -> message += "I don't know this code";
        }
        return message;
    }
}