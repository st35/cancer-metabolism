function r_ME2 = Malic_Enzyme_2(M_MAL, M_NAD, M_PYR, M_NADH, M_CO2, Rate)
    V_ME2_mf = 52.0*60.0; V_ME2_mr = 11.6*60.0;
    K_ME2_MAL = 0.4; K_ME2_NAD = 0.042;
    K_ME2_PYR = 45.0; K_ME2_NADH = 0.05; K_ME2_CO2 = 4.0;

    V_ME2_mf = V_ME2_mf*Rate;
    V_ME2_mr = V_ME2_mr*Rate;

    t0 = ((V_ME2_mf*M_MAL*M_NAD) / (K_ME2_MAL*K_ME2_NAD));
    t1 = ((V_ME2_mr*M_PYR*M_NADH*M_CO2) / (K_ME2_PYR*K_ME2_NADH*K_ME2_CO2));
    t2 = 1.0 + (M_MAL / K_ME2_MAL) + (M_NAD / K_ME2_NAD) + ((M_MAL*M_NAD) / (K_ME2_MAL*K_ME2_NAD));
    t3 = (M_PYR / K_ME2_PYR) + (M_NADH / K_ME2_NADH) + (M_CO2 / K_ME2_CO2);
    t4 = ((M_PYR*M_NADH) / (K_ME2_PYR*K_ME2_NADH)) + ((M_PYR*M_CO2) / (K_ME2_PYR*K_ME2_CO2)) + ((M_NADH*M_CO2) / (K_ME2_NADH*K_ME2_CO2));
    t5 = ((M_PYR*M_NADH*M_CO2) / (K_ME2_PYR*K_ME2_NADH*K_ME2_CO2));

    r_ME2 = (t0 - t1) / (t2 + t3 + t4 + t5);
end