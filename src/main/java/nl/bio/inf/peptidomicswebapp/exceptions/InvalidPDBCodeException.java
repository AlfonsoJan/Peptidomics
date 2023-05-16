package nl.bio.inf.peptidomicswebapp.exceptions;

/**
 * @Author: Seabarrel
 */
public class InvalidPDBCodeException extends Exception {

    public InvalidPDBCodeException() {
        super("Invalid PDB code.. needs to be 4 characters long, and exists in the database!");
    }

}
