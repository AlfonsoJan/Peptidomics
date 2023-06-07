package nl.bio.inf.peptidomicswebapp.exceptions;

/**
 * Custom exception for when the eigenvectors are not found for the api
 * @author Jan Alfonso Busker
 */
public class EigenVectorsNotFoundException extends RuntimeException {

    public EigenVectorsNotFoundException(Integer id) {
        super("Could not find eigenvectors with length:  " + id);
    }
}
