function r_MDH = Malate_Dehydrogenase(M_MAL, M_NAD, M_OAA, M_NADH, Rate)
    V_MDH_mf = 8.683*1e4*60.0; V_MDH_mr = 1.4*1e5*60.0;
    K_MDH_MAL = 987.0e-3; K_MDH_NAD = 538.0e-3;
    K_MDH_OAA = 1.58e-3; K_MDH_NADH = 0.59e-3;

    V_MDH_mf = V_MDH_mf*Rate;
    V_MDH_mr = V_MDH_mr*Rate;

    t0 = ((V_MDH_mf*M_MAL*M_NAD) / (K_MDH_MAL*K_MDH_NAD)) - ((V_MDH_mr*M_OAA*M_NADH) / (K_MDH_OAA*K_MDH_NADH));
    t1 = 1.0 + (M_MAL / K_MDH_MAL) + (M_NAD / K_MDH_NAD) + ((M_MAL*M_NAD) / (K_MDH_MAL*K_MDH_NAD));
    t2 = (M_OAA / K_MDH_OAA) + (M_NADH / K_MDH_NADH) + ((M_OAA*M_NADH) / (K_MDH_OAA*K_MDH_NADH));

    r_MDH = t0 / (t1 + t2);
end