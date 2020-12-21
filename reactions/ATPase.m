function r_ATP = ATPase(M_ATP, Psi, M_Pi, Rate)
    At = 4.16; K_app = 4.4e-6;
    Delta_G_transport = 1.2*96485.0*Psi; RT = 8.314*298.0;

    At = At*4.0;

    ATP_crit = (At) / (1.0 + (exp((-3.0*Delta_G_transport) / RT) / (K_app*M_Pi)));
    k_ATP = 131.9*3600.0; b = 4.0e-3;

    k_ATP = k_ATP*Rate;

    t0 = k_ATP;
    t1 = (2.0 / (1.0 + exp(b*(M_ATP - ATP_crit)))) - 1.0;

    r_ATP = t0*t1;
end