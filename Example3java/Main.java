//test example 3
public class Main{
    public static void main(String[] args) {
        // # of control variables
        final int noi = 7;

        // 7 control variables
        TPSAad x1 = new TPSAad();
        TPSAad x2 = new TPSAad();
        TPSAad x3 = new TPSAad();
        TPSAad x4 = new TPSAad();
        TPSAad x5 = new TPSAad();
        TPSAad x6 = new TPSAad();
        TPSAad x7 = new TPSAad();

        // 2D Twiss parameters (beta_x,alpha_x,beta_y,alpha_y)
        TPSAad[] twiss = new TPSAad[4];
        for (int i = 0; i < 4; i++) {
            twiss[i] = new TPSAad();
        }

        // Initialize each lattice variable
        x1.assign(0.2, 1);  // drift length
        x2.assign(0.1, 2);  // quad length
        x3.assign(29.6, 3); // quad strength (k)
        x4.assign(0.4, 4);  // drift length
        x5.assign(0.1, 5);  // quad length
        x6.assign(-29.6, 6);// quad strength (k)
        x7.assign(0.2, 7);  // drift length

        // initial Twiss parameters
        twiss[0].assign(1.0);
        twiss[1].assign(0.5);
        twiss[2].assign(1.0);
        twiss[3].assign(-0.5);

        // loop through nperd of the FODO lattice
        final int nperd = 1;
        for (int i = 0; i < nperd; i++) {
            // drift
            double xk = 0.0;
            quadmap(twiss, x1, x2, xk);
            // quad
            xk = x3.map[0];
            quadmap(twiss, x2, x3, xk);
            // drift
            xk = 0.0;
            quadmap(twiss, x4, x2, xk);
            // quad
            xk = x6.map[0];
            quadmap(twiss, x5, x6, xk);
            // drift
            xk = 0.0;
            quadmap(twiss, x7, x2, xk);
        }

        // Print results
        System.out.println("Beta_x and its derivatives w.r.t. 7 variables:");
        for (int i = 0; i < noi + 1; i++) {
            System.out.print(twiss[0].map[i] + " ");
        }
        System.out.println();

        System.out.println("Alpha_x and its derivatives w.r.t. 7 variables:");
        for (int i = 0; i < noi + 1; i++) {
            System.out.print(twiss[1].map[i] + " ");
        }
        System.out.println();

        System.out.println("Beta_y and its derivatives w.r.t. 7 variables:");
        for (int i = 0; i < noi + 1; i++) {
            System.out.print(twiss[2].map[i] + " ");
        }
        System.out.println();

        System.out.println("Alpha_y and its derivatives w.r.t. 7 variables:");
        for (int i = 0; i < noi + 1; i++) {
            System.out.print(twiss[3].map[i] + " ");
        }
        System.out.println();
    }

    private static void quadmap(TPSAad[] tpsaPtc, TPSAad tau, TPSAad outin, double xk) {
        TPSAad[] tpsaPtcTemp = new TPSAad[4];
        for (int i = 0; i < 4; i++) {
            tpsaPtcTemp[i] = new TPSAad();
        }

        TPSAad co = new TPSAad();
        TPSAad si = new TPSAad();
        TPSAad ch = new TPSAad();
        TPSAad sh = new TPSAad();
        TPSAad sqk = new TPSAad();
        
        // Calculate gammax and gammay
        TPSAad gammax = TPSAad.divide(
            TPSAad.add(1.0, TPSAad.pow(tpsaPtc[1], 2)),
            tpsaPtc[0]
        );
        
        TPSAad gammay = TPSAad.divide(
            TPSAad.add(1.0, TPSAad.pow(tpsaPtc[3], 2)),
            tpsaPtc[2]
        );

        if (xk == 0.0) {
            tpsaPtc[0] = TPSAad.subtract(
                tpsaPtc[0],
                TPSAad.subtract(
                    TPSAad.multiply(2.0, TPSAad.multiply(tpsaPtc[1], tau)),
                    TPSAad.multiply(tau, TPSAad.multiply(tau, gammax))
                )
            );
            
            tpsaPtc[2] = TPSAad.subtract(
                tpsaPtc[2],
                TPSAad.subtract(
                    TPSAad.multiply(2.0, TPSAad.multiply(tpsaPtc[3], tau)),
                    TPSAad.multiply(tau, TPSAad.multiply(tau, gammay))
                )
            );
        }
        else if (xk > 0.0) {
            sqk = TPSAad.sqrt(outin);
            co = TPSAad.cos(TPSAad.multiply(sqk, tau));
            si = TPSAad.sin(TPSAad.multiply(sqk, tau));
            ch = TPSAad.cosh(TPSAad.multiply(sqk, tau));
            sh = TPSAad.sinh(TPSAad.multiply(sqk, tau));

            // Calculate tpsaPtcTemp values
            tpsaPtcTemp[0] = calculateTpsaPtcTemp0(tpsaPtc, co, si, sqk, gammax);
            tpsaPtcTemp[1] = calculateTpsaPtcTemp1(tpsaPtc, co, si, sqk, gammax);
            tpsaPtcTemp[2] = calculateTpsaPtcTemp2(tpsaPtc, ch, sh, sqk, gammay);
            tpsaPtcTemp[3] = calculateTpsaPtcTemp3(tpsaPtc, ch, sh, sqk, gammay);

            System.arraycopy(tpsaPtcTemp, 0, tpsaPtc, 0, 4);
        }
        else if (xk < 0.0) {
            sqk = TPSAad.sqrt(TPSAad.multiply(-1.0, outin));
            co = TPSAad.cos(TPSAad.multiply(sqk, tau));
            si = TPSAad.sin(TPSAad.multiply(sqk, tau));
            ch = TPSAad.cosh(TPSAad.multiply(sqk, tau));
            sh = TPSAad.sinh(TPSAad.multiply(sqk, tau));

            // Calculate tpsaPtcTemp values for negative xk
            tpsaPtcTemp[0] = calculateNegativeTpsaPtcTemp0(tpsaPtc, ch, sh, sqk, gammax);
            tpsaPtcTemp[1] = calculateNegativeTpsaPtcTemp1(tpsaPtc, ch, sh, sqk, gammax);
            tpsaPtcTemp[2] = calculateNegativeTpsaPtcTemp2(tpsaPtc, co, si, sqk, gammay);
            tpsaPtcTemp[3] = calculateNegativeTpsaPtcTemp3(tpsaPtc, co, si, sqk, gammay);

            System.arraycopy(tpsaPtcTemp, 0, tpsaPtc, 0, 4);
        }
    }

    // Helper methods for calculating tpsaPtcTemp values
    // ... (include the helper methods from the previous response)
    private static TPSAad calculateTpsaPtcTemp0(TPSAad[] tpsaPtc, 
            TPSAad co, TPSAad si, TPSAad sqk, TPSAad gammax) {
        
        TPSAad term1 = TPSAad.multiply(tpsaPtc[0], TPSAad.multiply(co, co));
        
        TPSAad term2 = TPSAad.multiply(2.0, 
                       TPSAad.multiply(tpsaPtc[1], 
                       TPSAad.multiply(co, 
                       TPSAad.divide(si, sqk))));
        
        TPSAad siOverSqk = TPSAad.divide(si, sqk);
        TPSAad term3 = TPSAad.multiply(
                       TPSAad.multiply(siOverSqk, siOverSqk), 
                       gammax);
        
        return TPSAad.subtract(term1, TPSAad.subtract(term2, term3));
    }

    private static TPSAad calculateTpsaPtcTemp1(TPSAad[] tpsaPtc, 
            TPSAad co, TPSAad si, TPSAad sqk, TPSAad gammax) {
        
        TPSAad term1 = TPSAad.multiply(tpsaPtc[0], 
                       TPSAad.multiply(TPSAad.multiply(co, si), sqk));
        
        TPSAad coSquared = TPSAad.multiply(co, co);
        TPSAad siSquared = TPSAad.multiply(si, si);
        TPSAad term2 = TPSAad.multiply(tpsaPtc[1], 
                       TPSAad.subtract(coSquared, siSquared));
        
        TPSAad term3 = TPSAad.multiply(TPSAad.multiply(co, si), 
                       TPSAad.divide(gammax, sqk));
        
        return TPSAad.add(term1, TPSAad.subtract(term2, term3));
    }

    private static TPSAad calculateTpsaPtcTemp2(TPSAad[] tpsaPtc, 
            TPSAad ch, TPSAad sh, TPSAad sqk, TPSAad gammay) {
        
        TPSAad term1 = TPSAad.multiply(tpsaPtc[2], 
                       TPSAad.multiply(ch, ch));
        
        TPSAad term2 = TPSAad.multiply(2.0, 
                       TPSAad.multiply(tpsaPtc[3], 
                       TPSAad.multiply(ch, 
                       TPSAad.divide(sh, sqk))));
        
        TPSAad shOverSqk = TPSAad.divide(sh, sqk);
        TPSAad term3 = TPSAad.multiply(
                       TPSAad.multiply(shOverSqk, shOverSqk), 
                       gammay);
        
        return TPSAad.subtract(term1, TPSAad.subtract(term2, term3));
    }

    private static TPSAad calculateTpsaPtcTemp3(TPSAad[] tpsaPtc, 
            TPSAad ch, TPSAad sh, TPSAad sqk, TPSAad gammay) {
        
        TPSAad term1 = TPSAad.multiply(TPSAad.multiply(-1.0, tpsaPtc[2]), 
                       TPSAad.multiply(TPSAad.multiply(ch, sh), sqk));
        
        TPSAad chSquared = TPSAad.multiply(ch, ch);
        TPSAad shSquared = TPSAad.multiply(sh, sh);
        TPSAad term2 = TPSAad.multiply(tpsaPtc[3], 
                       TPSAad.add(chSquared, shSquared));
        
        TPSAad term3 = TPSAad.multiply(TPSAad.multiply(ch, sh), 
                       TPSAad.divide(gammay, sqk));
        
        return TPSAad.add(term1, TPSAad.subtract(term2, term3));
    }

   private static TPSAad calculateNegativeTpsaPtcTemp0(TPSAad[] tpsaPtc,
            TPSAad ch, TPSAad sh, TPSAad sqk, TPSAad gammax) {

        TPSAad term1 = TPSAad.multiply(tpsaPtc[0],
                       TPSAad.multiply(ch, ch));

        TPSAad term2 = TPSAad.multiply(2.0,
                       TPSAad.multiply(tpsaPtc[1],
                       TPSAad.multiply(ch,
                       TPSAad.divide(sh, sqk))));

        TPSAad shOverSqk = TPSAad.divide(sh, sqk);
        TPSAad term3 = TPSAad.multiply(
                       TPSAad.multiply(shOverSqk, shOverSqk),
                       gammax);

        return TPSAad.subtract(term1, TPSAad.subtract(term2, term3));
    }

    private static TPSAad calculateNegativeTpsaPtcTemp1(TPSAad[] tpsaPtc,
            TPSAad ch, TPSAad sh, TPSAad sqk, TPSAad gammax) {

        TPSAad term1 = TPSAad.multiply(TPSAad.multiply(-1.0, tpsaPtc[0]),
                       TPSAad.multiply(TPSAad.multiply(ch, sh), sqk));

        TPSAad chSquared = TPSAad.multiply(ch, ch);
        TPSAad shSquared = TPSAad.multiply(sh, sh);
        TPSAad term2 = TPSAad.multiply(tpsaPtc[1],
                       TPSAad.add(chSquared, shSquared));

        TPSAad term3 = TPSAad.multiply(TPSAad.multiply(ch, sh),
                       TPSAad.divide(gammax, sqk));

        return TPSAad.add(term1, TPSAad.subtract(term2, term3));
    }

    private static TPSAad calculateNegativeTpsaPtcTemp2(TPSAad[] tpsaPtc,
            TPSAad co, TPSAad si, TPSAad sqk, TPSAad gammay) {

        TPSAad term1 = TPSAad.multiply(tpsaPtc[2], TPSAad.multiply(co, co));

        TPSAad term2 = TPSAad.multiply(2.0,
                       TPSAad.multiply(tpsaPtc[3],
                       TPSAad.multiply(co,
                       TPSAad.divide(si, sqk))));

        TPSAad siOverSqk = TPSAad.divide(si, sqk);
        TPSAad term3 = TPSAad.multiply(
                       TPSAad.multiply(siOverSqk, siOverSqk),
                       gammay);

        return TPSAad.subtract(term1, TPSAad.subtract(term2, term3));
    }

    private static TPSAad calculateNegativeTpsaPtcTemp3(TPSAad[] tpsaPtc,
            TPSAad co, TPSAad si, TPSAad sqk, TPSAad gammay) {

        TPSAad term1 = TPSAad.multiply(tpsaPtc[2],
                       TPSAad.multiply(TPSAad.multiply(co, si), sqk));

        TPSAad coSquared = TPSAad.multiply(co, co);
        TPSAad siSquared = TPSAad.multiply(si, si);
        TPSAad term2 = TPSAad.multiply(tpsaPtc[3],
                       TPSAad.subtract(coSquared, siSquared));

        TPSAad term3 = TPSAad.multiply(TPSAad.multiply(co, si),
                       TPSAad.divide(gammay, sqk));

        return TPSAad.add(term1, TPSAad.subtract(term2, term3));
    }

}
