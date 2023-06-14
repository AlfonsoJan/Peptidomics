package nl.bio.inf.peptidomicswebapp.config;

import jakarta.servlet.http.HttpSessionEvent;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import org.springframework.security.web.session.HttpSessionEventPublisher;
import org.springframework.stereotype.Component;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.logging.Logger;

/**
 *  This class handles the creation and deletion of the session.
 * @author Jan Alfonso
 */
@Component
public class SessionDestroyer extends HttpSessionEventPublisher {

    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());

    /**
     * Creates the session
     * @param event HttpSessionEvent passed in by the container
     */
    @Override
    public void sessionCreated(HttpSessionEvent event) {
        LOGGER.info("Created session: " + event.getSession().getId());
        super.sessionCreated(event);
    }

    /**
     * Destroys the session and deletes the files in the session
     * @param event The HttpSessionEvent pass in by the container
     */
    @Override
    public void sessionDestroyed(HttpSessionEvent event) {
        LOGGER.info("Session destroyed: " + event.getSession().getId());
        deleteTempFiles(event);
        super.sessionDestroyed(event);
    }

    /**
     * Deletes the temporary files from the system!
     * @param event session event
     * @throws RuntimeException when files can not be accessed
     */
    public void deleteTempFiles(HttpSessionEvent event) {
        String tempLocation = String.valueOf(event.getSession().getAttribute("tempLocation"));
        String jsonFileLocation = String.valueOf(event.getSession().getAttribute("jsonFile"));
        try {
            Files.delete(Path.of(tempLocation));
            Files.delete(Path.of(jsonFileLocation));
            LOGGER.info("Deleted session files of: " + event.getSession().getId());
        } catch (IOException e) {
            LOGGER.warning("Could not delete the files of session: " + event.getSession().getId());
        }
    }
}
