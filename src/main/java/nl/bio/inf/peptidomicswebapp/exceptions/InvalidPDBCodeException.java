package nl.bio.inf.peptidomicswebapp.exceptions;

/**
 * Custom exception for when the given pdb code is not a valid one
 * @author Wouter Zeevat
 */
public class InvalidPDBCodeException extends Exception {

    public InvalidPDBCodeException() {
        super("Invalid PDB code.. needs to be 4 characters long, and exists in the database!");
    }

}
