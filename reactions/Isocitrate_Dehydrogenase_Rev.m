function r_IDH_Rev = Isocitrate_Dehydrogenase_Rev(M_NADPH, M_AKG, M_CO2, M_CIT, M_NADP, Rate)
    V_IDHcP_mf = 4.25e-1*1e3*3600.0; V_IDHcP_mr = 0.0;

    K_IDHcP_NADPH = 35.0e-3; K_IDHcP_AKG = 1.4; K_IDHcP_CO2 = 12.6;
    K_IDHcP_CIT = 5.0e-4; K_IDHcP_NADP = 1.0e-4;

    V_IDHcP_mf = V_IDHcP_mf*Rate;

    t0 = ((V_IDHcP_mf*M_NADPH*M_AKG*M_CO2) / (K_IDHcP_NADPH*K_IDHcP_AKG*K_IDHcP_CO2));
    t1 = ((V_IDHcP_mr*M_CIT*M_NADP) / (K_IDHcP_CIT*K_IDHcP_NADP));
    t2 = 1.0 + (M_NADPH / K_IDHcP_NADPH) + (M_AKG / K_IDHcP_AKG) + (M_CO2 / K_IDHcP_CO2);
    t3 = ((M_NADPH*M_AKG) / (K_IDHcP_NADPH*K_IDHcP_AKG)) + ((M_NADPH*M_CO2) / (K_IDHcP_NADPH*K_IDHcP_CO2)) + ((M_AKG*M_CO2) / (K_IDHcP_AKG*K_IDHcP_CO2));
    t4 = ((M_NADPH*M_AKG*M_CO2) / (K_IDHcP_NADPH*K_IDHcP_AKG*K_IDHcP_CO2));
    t5 = (M_CIT / K_IDHcP_CIT) + (M_NADP / K_IDHcP_NADP);
    t6 = ((M_CIT*M_NADP) / (K_IDHcP_CIT*K_IDHcP_NADP));

    r_IDH_Rev = (t0 - t1) / (t2 + t3 + t4 + t5 + t6);
end