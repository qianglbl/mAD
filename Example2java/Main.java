public class Main {
    private static final int NOI = 7;  // # of control variables

    public static void main(String[] args) {
        // 7 control variables
        TPSAad x1 = new TPSAad();
        TPSAad x2 = new TPSAad();
        TPSAad x3 = new TPSAad();
        TPSAad x4 = new TPSAad();
        TPSAad x5 = new TPSAad();
        TPSAad x6 = new TPSAad();
        TPSAad x7 = new TPSAad();

        // particle coordinates (x,px,y,py)
        TPSAad[] pt = new TPSAad[4];
        for (int i = 0; i < 4; i++) {
            pt[i] = new TPSAad();
        }

        // Initialize each lattice variable
        x1.assign(0.2, 1);    // drift length
        x2.assign(0.1, 2);    // quad length
        x3.assign(29.6, 3);   // quad strength (k)
        x4.assign(0.4, 4);    // drift length
        x5.assign(0.1, 5);    // quad length
        x6.assign(-29.6, 6);  // quad strength (k)
        x7.assign(0.2, 7);    // drift length

        // initial particle coordinates
        pt[0].assign(1.0e-3, 0);  // x
        pt[1].assign(1.0e-3, 0);  // px
        pt[2].assign(-1.0e-3, 0); // y
        pt[3].assign(-1.0e-3, 0); // py

        // loop through nperd of the FODO lattice
        final int nperd = 1;
        for (int i = 0; i < nperd; i++) {
            // drift
            double xk = 0.0;
            quadmap(pt, x1, x2, xk);
            // quad
            xk = x3.map[0];
            quadmap(pt, x2, x3, xk);
            // drift
            xk = 0.0;
            quadmap(pt, x4, x2, xk);
            // quad
            xk = x6.map[0];
            quadmap(pt, x5, x6, xk);
            // drift
            xk = 0.0;
            quadmap(pt, x7, x2, xk);
        }

        // Print results
        System.out.println("X coordinate and its derivatives w.r.t. 7 variables:");
        for (int i = 0; i < NOI + 1; i++) {
            System.out.print(pt[0].map[i] + " ");
        }
        System.out.println();

        System.out.println("Px coordinate and its derivatives w.r.t. 7 variables:");
        for (int i = 0; i < NOI + 1; i++) {
            System.out.print(pt[1].map[i] + " ");
        }
        System.out.println();

        System.out.println("Y coordinate and its derivatives w.r.t. 7 variables:");
        for (int i = 0; i < NOI + 1; i++) {
            System.out.print(pt[2].map[i] + " ");
        }
        System.out.println();

        System.out.println("Py coordinate and its derivatives w.r.t. 7 variables:");
        for (int i = 0; i < NOI + 1; i++) {
            System.out.print(pt[3].map[i] + " ");
        }
        System.out.println();
    }

    private static void quadmap(TPSAad[] tpsaPtc, TPSAad tau, TPSAad outin, double xk) {
        TPSAad[] tpsaPtcTemp = new TPSAad[4];
        for (int i = 0; i < 4; i++) {
            tpsaPtcTemp[i] = new TPSAad();
        }

        if (xk == 0.0) {
            // Drift section
            tpsaPtc[0] = tpsaPtc[0].add(tpsaPtc[1].multiply(tau));
            tpsaPtc[2] = tpsaPtc[2].add(tpsaPtc[3].multiply(tau));
        }
        else if (xk > 0.0) {
            TPSAad sqk = TPSAad.sqrt(outin);
            TPSAad tauSqk = sqk.multiply(tau);
            TPSAad co = TPSAad.cos(tauSqk);
            TPSAad si = TPSAad.sin(tauSqk);
            TPSAad ch = TPSAad.cosh(tauSqk);
            TPSAad sh = TPSAad.sinh(tauSqk);

            tpsaPtcTemp[0] = tpsaPtc[0].multiply(co).add(
                            tpsaPtc[1].multiply(si).divide(sqk));
            tpsaPtcTemp[1] = tpsaPtc[0].multiply(si).multiply(sqk).multiply(new TPSAad(-1.0)).add(
                            tpsaPtc[1].multiply(co));
            tpsaPtcTemp[2] = tpsaPtc[2].multiply(ch).add(
                            tpsaPtc[3].multiply(sh).divide(sqk));
            tpsaPtcTemp[3] = tpsaPtc[2].multiply(sh).multiply(sqk).add(
                            tpsaPtc[3].multiply(ch));

            System.arraycopy(tpsaPtcTemp, 0, tpsaPtc, 0, 4);
        }
        else if (xk < 0.0) {
            TPSAad sqk = TPSAad.sqrt(outin.multiply(new TPSAad(-1.0)));
            TPSAad tauSqk = sqk.multiply(tau);
            TPSAad co = TPSAad.cos(tauSqk);
            TPSAad si = TPSAad.sin(tauSqk);
            TPSAad ch = TPSAad.cosh(tauSqk);
            TPSAad sh = TPSAad.sinh(tauSqk);

            tpsaPtcTemp[0] = tpsaPtc[0].multiply(ch).add(
                            tpsaPtc[1].multiply(sh).divide(sqk));
            tpsaPtcTemp[1] = tpsaPtc[0].multiply(sh).multiply(sqk).add(
                            tpsaPtc[1].multiply(ch));
            tpsaPtcTemp[2] = tpsaPtc[2].multiply(co).add(
                            tpsaPtc[3].multiply(si).divide(sqk));
            tpsaPtcTemp[3] = tpsaPtc[2].multiply(si).multiply(sqk).multiply(new TPSAad(-1.0)).add(
                            tpsaPtc[3].multiply(co));

            System.arraycopy(tpsaPtcTemp, 0, tpsaPtc, 0, 4);
        }
    }
}
