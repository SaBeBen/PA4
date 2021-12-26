package ode;

/**
 * Das Einschrittverfahren von Heun
 *
 * @author braeckle
 *
 */
public class Heun implements Einschrittverfahren {

    @Override
    /**
     * {@inheritDoc}
     * Nutzen Sie dabei geschickt den Expliziten Euler.
     */
    public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
        // TODO: diese Methode ist zu implementieren
        double[] result = new double[y_k.length];
        ExpliziterEuler euler = new ExpliziterEuler();
        double[] ausgewertet = ode.auswerten(t, y_k);
        double[] ausgewertetEuler = ode.auswerten(t + delta_t, euler.nextStep(y_k, t, delta_t, ode));
        for (int i = 0; i < result.length; i++) {
            result[i] = y_k[i] + (delta_t/2)*(ausgewertet[i] + ausgewertetEuler[i]);
        }
        return result;
    }

}
