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
     * Destroys the session
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
     * @param event
     */
    private void deleteTempFiles(HttpSessionEvent event) {
        String tempLocationCompare = String.valueOf(event.getSession().getAttribute("tempLocationCompare"));
        String tempLocation = String.valueOf(event.getSession().getAttribute("tempLocation"));
        try {
            Files.delete(Path.of(tempLocationCompare));
            Files.delete(Path.of(tempLocation));
            LOGGER.info("Deleted session files of: " + event.getSession().getId());
        } catch (IOException e) {
            LOGGER.warning("Could not delete the files!!");
            throw new RuntimeException(e);
        }
    }
}
