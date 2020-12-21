function r_FUM = Fumarase(M_FUM, M_MAL, M_CIT, M_ATP, M_ADP, M_GTP, M_GDP, Rate)
    V_FUM_mf = 2.712*1e5*60.0; V_FUM_mr = 2.736*1e5*60.0;
    K_FUM_FUM = 44.7e-3; K_FUM_MAL = 197.7e-3;
    K_FUM_CITi = 3500.0e-3; K_FUM_ATPi = 40.0e-3; K_FUM_ADPi = 400.0e-3; K_FUM_GTPi = 80.0e-3; K_FUM_GDPi = 330.0e-3;

    V_FUM_mf = V_FUM_mf*Rate;
    V_FUM_mr = V_FUM_mr*Rate;

    t0 = ((V_FUM_mf*M_FUM) / K_FUM_FUM) - ((V_FUM_mr*M_MAL) / K_FUM_MAL);
    t1 = 1.0 + (M_FUM / K_FUM_FUM) + (M_MAL / K_FUM_MAL);
    alpha_i = (M_CIT / K_FUM_CITi) + (M_ATP / K_FUM_ATPi) + (M_ADP / K_FUM_ADPi) + (M_GTP / K_FUM_GTPi) + (M_GDP / K_FUM_GDPi);

    r_FUM = t0 / (alpha_i + t1);
end