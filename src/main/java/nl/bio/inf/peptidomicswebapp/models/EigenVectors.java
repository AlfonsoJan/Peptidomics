package nl.bio.inf.peptidomicswebapp.models;

import java.util.ArrayList;
import java.util.List;


public class EigenVectors {

    public int length;
    public List<String> x;
    public List<String> y;
    public List<String> z;

    public void addLine(String csvLine) {
        String[] elements = csvLine.split(",");
        x.add(elements[0]);
        y.add(elements[1]);
        z.add(elements[2]);
    }

    public EigenVectors(Integer length) {
        if (length == null) throw new NullPointerException();
        this.length = length;
        this.x = new ArrayList<>();
        this.y = new ArrayList<>();
        this.z = new ArrayList<>();
    }

    @Override
    public String toString() {
        return "EigenVectors{" +
                "x=" + x +
                ", y=" + y +
                ", z=" + z +
                ", length=" + length +
                '}';
    }

    public List<String> getX() {
        return x;
    }

    public List<String> getY() {
        return y;
    }

    public List<String> getZ() {
        return z;
    }

    public int getLength() {
        return length;
    }
}
