package nl.bio.inf.peptidomicswebapp.models;

import nl.bio.inf.peptidomicswebapp.exceptions.EigenVectorsNotFoundException;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.junit.jupiter.api.Assertions.*;

public class EigenVectorsTest {

    @Test
    void createEigenVector() {
        EigenVectors vector = new EigenVectors(1);
        assertEquals(1, vector.getLength());
    }

    @Test
    void createEigenVectorNull() {
        assertThrows(NullPointerException.class, () -> new EigenVectors(null));
    }

    @Test
    void addLine() {
        EigenVectors vector = new EigenVectors(1);
        vector.addLine("1,0,0");
        vector.addLine("1,0,0");
        assertEquals(vector.getX().size(), 2);
    }

    @ParameterizedTest
    @ValueSource(ints = {1,2,3,41,5,2,52})
    void testAddlineMultiple(int integer) {
        EigenVectors vector = new EigenVectors(1);
        vector.addLine(String.valueOf(integer) + ",0,0");
        assertEquals(vector.getX().get(0), String.valueOf(integer));
    }

}
