//test example 1
public class Main {
    public static void main(String[] args) {
        TPSAad x1 = new TPSAad();
        TPSAad x2 = new TPSAad();
        TPSAad x3 = new TPSAad();
        TPSAad x4 = new TPSAad();
        TPSAad x5 = new TPSAad();
        TPSAad func;

        x1.assign(3.0, 1);
        x2.assign(0.1, 2);
        x3.assign(0.1, 3);
        x4.assign(1.0, 4);
        x5.assign(1.0, 5);

        // test function
        // func = (2*cos(x1/x2)+x1/x2+exp(x2))*x3+2*x4+sinh(x5);
        TPSAad term1 = TPSAad.cos(x1.divide(x2)).multiply(new TPSAad(2.0));
        TPSAad term2 = x1.divide(x2);
        TPSAad term3 = TPSAad.exp(x2);
        TPSAad term4 = x4.multiply(new TPSAad(2.0));
        TPSAad term5 = TPSAad.sinh(x5);

        func = term1.add(term2).add(term3).multiply(x3).add(term4).add(term5);

        System.out.println("func. value and its derivatives w.r.t. 5 variables");
        System.out.println(func.map[0] + " " + func.map[1] + " " + func.map[2] + " " + 
                         func.map[3] + " " + func.map[4] + " " + func.map[5]);
    }
}
