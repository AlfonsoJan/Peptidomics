package nl.bio.inf.peptidomicswebapp.exceptions;

/**
 * @Author: Seabarrel
 */
public class TooLargeNumberException extends Exception {

    public TooLargeNumberException() {
        super("This number is too large!");
    }

}
