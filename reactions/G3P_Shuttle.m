function r_G3P_Shuttle = G3P_Shuttle(C_G3P, Psi, O2, Rate)
    k_resp = 2.5*3600.0; K = 0.002e-3; a = 100.0; Psi_m = 0.15;
    k_O2 = 0.1;

    K = 0.01;

    k_resp = k_resp*Rate;

    t1 = 1.0 / (1.0 + exp(a*(Psi - Psi_m)));
    t2 = (O2) / (k_O2 + O2);

    K_G3P_Shuttle_G3P = 1.2;
    t3 = (k_resp*C_G3P) / (K_G3P_Shuttle_G3P + C_G3P);

    r_G3P_Shuttle = t1*t2*t3;
end