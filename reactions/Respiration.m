function r_Resp = Respiration(M_NADH, Psi, O2, Rate)
    k_resp = 2.5*3600.0; K = 0.002e-3; a = 100.0; Psi_m = 0.15;
    k_O2 = 0.1;

    K = 0.01;

    k_resp = k_resp*Rate;

    t0 = (k_resp*M_NADH) / (K + M_NADH);
    t1 = 1.0 / (1.0 + exp(a*(Psi - Psi_m)));
    t2 = (O2) / (k_O2 + O2);

    r_Resp = t0*t1*t2;
end