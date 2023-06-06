package nl.bio.inf.peptidomicswebapp.exceptions;

public class EigenVectorsNotFoundException extends RuntimeException {

    public EigenVectorsNotFoundException(Integer id) {
        super("Could not find eigenvectors with length:  " + id);
    }
}
