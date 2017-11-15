package net.yuri.klymovych.nonlinearfitter.classes;

public class FitHost implements LMhost{
    //--------constants for LM---------------
    private final double DELTAP = 1e-6; // parm step
    private final double BIGVAL = 9e99; // fault flag

    //---constants for fitting: pixels & parameters-----

    private double data[][] =   // {x, y} data pairs [row][col]
            {{0.00,	0.6793},
                    {0.03,	0.6787},
                    {0.06,	0.6768},
                    {0.09,	0.6736},
                    {0.12,	0.6691},
                    {0.15,	0.6634},
                    {0.18,	0.6565},
                    {0.21,	0.6482},
                    {0.24,	0.6388},
                    {0.27,	0.6280},
                    {0.30,	0.6161},
                    {0.33,	0.6030},
                    {0.36,	0.5887},
                    {0.39,	0.5733},
                    {0.42,	0.5568},
                    {0.45,	0.5394},
                    {0.48,	0.5210},
                    {0.51,	0.5019},
                    {0.54,	0.4820},
                    {0.57,	0.4614},
                    {0.60,	0.4404},
                    {0.63,	0.4191},
                    {0.66,	0.3975},
                    {0.69,	0.3758},
                    {0.72,	0.3542},
                    {0.75,	0.3328},
                    {0.78,	0.3117},
                    {0.81,	0.2910},
                    {0.84,	0.2710},
                    {0.87, 0.2515},
                    {0.90,	0.2328},
                    {0.93,	0.2149},
                    {0.96,	0.1979},
                    {0.99,	0.1818},
                    {1.02,	0.1667},
                    {1.05,	0.1524},
                    {1.08,	0.1392},
                    {1.11,	0.1268},
                    {1.14,	0.1154},
                    {1.17,	0.1048},
                    {1.20,	0.0951},
                    {1.23,	0.0862},
                    {1.26,	0.0781},
                    {1.29,	0.0706},
                    {1.32,	0.0639},
                    {1.35,	0.0577},
                    {1.38,	0.0521},
                    {1.41,	0.0471},
                    {1.44,	0.0425},
                    {1.47,	0.0384},
                    {1.50,	0.0347},
                    {1.53,	0.0313},
                    {1.56,	0.0283},
                    {1.59,	0.0255},
                    {1.62,	0.0231},
                    {1.65,	0.0209},
                    {1.68,	0.0189},
                    {1.71,	0.0171},
                    {1.74,	0.0155},
                    {1.77,	0.0140},
                    {1.80,	0.0127},
                    {1.83,	0.0115},
                    {1.86,	0.0105},
                    {1.89,	0.0095},
                    {1.92,	0.0087},
                    {1.95,	0.0079},
                    {1.98,	0.0072}};

    private int    NPTS     = data.length;
    private double WEIGHT   = 1.0;
    private double parms[]  = {1.0, 1.0, 1.0, 1.0}; // initial guess
    private int    NPARMS   = parms.length;

    private double   resid[]  = new double[NPTS];
    private double jac[][] = new double[NPTS][NPARMS];

    void doFit() {
        listParms("Start Parms", parms);
        LM myLM = new LM(this, NPARMS, NPTS);
        listParms("Fitted Parms", parms);
    }

   /* public FitHost()
    {
        listParms("Start Parms", parms);
        LM myLM = new LM(this, NPARMS, NPTS);
        listParms("Fitted Parms", parms);
    }*/

    //------------mathematical helpers for FitHost--------------


    private void listParms(String title, double p[])
    {
        System.out.print(title + "----");
        for (int i=0; i<NPARMS; i++)
            System.out.print(String.format("%12.6f", p[i]));
        System.out.println("");
    }


    private double func(int i, double data[][], double p[])
    // inverse even polynomial.
    // Called only by dComputeResiduals().
    {
        double x = data[i][0];
        double x2 = x*x;
        double x4 = x2*x2;
        double x6 = x4*x2;
        // double x8 = x6*x2;
        double denom = 1 + p[1]*x2 + p[2]*x4 + p[3]*x6;  // + p[4]*x8;
        return p[0]/denom;
    }


    private double dComputeResiduals()
    // Called by mandated dNudge().
    // Uses current parms[] vector;
    // Evaluates resid[i] = (model[i] - data[i])*WEIGHT.
    // Returns sum-of-squares.
    {
        double sumsq = 0.0;
        for (int i=0; i<NPTS; i++)
        {
            double y = data[i][1];  // row i, col 1
            resid[i] = (func(i,data,parms) - y) * WEIGHT;
            sumsq += resid[i]*resid[i];
        }
        return sumsq;
    }

    //------the four mandated interface methods------------

    public double dNudge(double dp[])
    // Allows LM to modify parms[] and reevaluate its fit.
    // Returns sum-of-squares for nudged params.
    // This is the only place that parms[] are modified.
    // If NADJ<NPARMS, this is the place for your LUT.
    {
        for (int j=0; j<NPARMS; j++)
            parms[j] += dp[j];
        return dComputeResiduals();
    }

    public boolean bBuildJacobian()
    // Allows LM to compute a new Jacobian.
    // Uses current parms[] and two-sided finite difference.
    // If current parms[] is bad, returns false.
    {
        double delta[] = new double[NPARMS];
        double FACTOR = 0.5 / DELTAP;
        double d;

        for (int j=0; j<NPARMS; j++)
        {
            for (int k=0; k<NPARMS; k++)
                delta[k] = (k==j) ? DELTAP : 0.0;

            d = dNudge(delta); // resid at pplus
            if (d==BIGVAL)
            {
                System.out.println("Bad dBuildJacobian() exit 2");
                return false;
            }
            for (int i=0; i<NPTS; i++)
                jac[i][j] = dGetResid(i);

            for (int k=0; k<NPARMS; k++)
                delta[k] = (k==j) ? -2*DELTAP : 0.0;

            d = dNudge(delta); // resid at pminus
            if (d==BIGVAL)
            {
                System.out.println("Bad dBuildJacobian() exit 3");
                return false;
            }

            for (int i=0; i<NPTS; i++)
                jac[i][j] -= dGetResid(i);  // fetches resid[]

            for (int i=0; i<NPTS; i++)
                jac[i][j] *= FACTOR;

            for (int k=0; k<NPARMS; k++)
                delta[k] = (k==j) ? DELTAP : 0.0;

            d = dNudge(delta);
            if (d==BIGVAL)
            {
                System.out.println("Bad dBuildJacobian() exit 4");
                return false;
            }
        }
        return true;
    }

    public double dGetResid(int i)
    // Allows LM to see one element of the resid[] vector.
    {
        return resid[i];
    }

    public double dGetJac(int i, int j)
    // Allows LM to see one element of the Jacobian matrix.
    {
        return jac[i][j];
    }
}
