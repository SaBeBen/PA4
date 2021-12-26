package ode;

/**
 * Der klassische Runge-Kutta der Ordnung 4
 *
 * @author braeckle
 */
public class RungeKutta4 implements Einschrittverfahren {

    @Override
    /**
     * {@inheritDoc}
     * Bei der Umsetzung koennen die Methoden addVectors und multScalar benutzt werden.
     */
    public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {

        double[] k1 = new double[y_k.length];
        double[] k2 = new double[y_k.length];
        double[] k3 = new double[y_k.length];
        double[] k4 = new double[y_k.length];
        double[] temp = new double[y_k.length];

        double[] k1Ausgewertet = ode.auswerten(t, y_k);
        for (int i = 0; i < y_k.length; i++) {
            k1[i] = delta_t * k1Ausgewertet[i];
            temp[i] = k1[i] / 2;
        }
        // k2 berechnen
        for (int i = 0; i < y_k.length; i++) {
            temp[i] = temp[i] + y_k[i];
        }
        double[] k2Ausgewertet = ode.auswerten(t + delta_t / 2, temp);
        for (int i = 0; i < y_k.length; i++) {
            k2[i] = delta_t * k2Ausgewertet[i];
            temp[i] = k2[i] / 2;
        }
        // k3 berechnen
        for (int i = 0; i < y_k.length; i++) {
            temp[i] = temp[i] + y_k[i];
        }
        double[] k3Ausgewertet = ode.auswerten(t + delta_t / 2, temp);
        for (int i = 0; i < y_k.length; i++) {
            k3[i] = delta_t * k3Ausgewertet[i];
            temp[i] = k3[i];
        }
        // k4 berechnen
        for (int i = 0; i < y_k.length; i++) {
            temp[i] = temp[i] + y_k[i];
        }
        double[] k4Ausgewertet = ode.auswerten(t + delta_t, temp);
        for (int i = 0; i < y_k.length; i++) {
            k4[i] = delta_t * k4Ausgewertet[i];
        }

        double[] result = new double[y_k.length];
        for (int i = 0; i < y_k.length; i++) {
            result[i] = y_k[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
        }

        return result;
    }

    /**
     * addiert die zwei Vektoren a und b
     */
    private double[] addVectors(double[] a, double[] b) {
        double[] erg = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            erg[i] = a[i] + b[i];
        }
        return erg;
    }

    /**
     * multipliziert den Skalar scalar auf den Vektor a
     */
    private double[] multScalar(double[] a, double scalar) {
        double[] erg = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            erg[i] = scalar * a[i];
        }
        return erg;
    }

}
