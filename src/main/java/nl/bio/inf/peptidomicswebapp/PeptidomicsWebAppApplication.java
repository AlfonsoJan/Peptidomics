package nl.bio.inf.peptidomicswebapp;

import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Logger;

import nl.bio.inf.peptidomicswebapp.config.HtmlLogFormatter;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.autoconfigure.security.servlet.UserDetailsServiceAutoConfiguration;

/**
 * The main to run the spring boot project
 * @author Wouter Zeevat
 * @author Jan Alfonso Busker
 */

@SpringBootApplication(exclude= {UserDetailsServiceAutoConfiguration.class})
public class PeptidomicsWebAppApplication {
    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());

    public static void main(String[] args) throws IOException {
        SpringApplication.run(PeptidomicsWebAppApplication.class, args);

        // Set the file handler for the logger-

        Handler[] handlers = LOGGER.getHandlers();
        for(Handler handler : handlers) {
            LOGGER.removeHandler(handler);
        }
        FileHandler fileHandler = new FileHandler("PeptidomicsLog.html", false);
        fileHandler.setFormatter(new HtmlLogFormatter());
        LOGGER.addHandler(fileHandler);
        LOGGER.info("--- Application Started! ---");
    }
}
