package nl.bio.inf.peptidomicswebapp.exceptions;

/**
 * Custom exception for when the oligo length is larger than set in the application properties
 * @author Wouter Zeevat
 */
public class TooLargeNumberException extends Exception {

    public TooLargeNumberException() {
        super("This number is too large!");
    }

}
