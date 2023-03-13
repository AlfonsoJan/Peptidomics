package nl.bio.inf.peptidomicswebapp.service;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.PropertySource;
import org.springframework.core.env.Environment;
import org.springframework.stereotype.Service;

@Service
@PropertySource("classpath:httpErrorCodes.properties")
public class ErrorService {

    @Autowired
    private Environment env;

    public String generateErrorMessage(final int errorCode) {
        String message = "There was a code " + errorCode + " error: ";
        switch (errorCode) {
            case 302:
                message += env.getProperty("302");
                break;
            case 304:
                message += env.getProperty("304");
                break;
            case 400:
                message += env.getProperty("400");
                break;
            case 401:
                message += env.getProperty("401");
                break;
            case 404:
                message += env.getProperty("404");
                break;
            case 405:
                message += env.getProperty("405");
                break;
            case 410:
                message += env.getProperty("410");
                break;
            case 500:
                message += env.getProperty("500");
                break;
            case 502:
                message += env.getProperty("502");
                break;
            default:
                message += "I don't know this code";
        }
        return message;
    }
}
