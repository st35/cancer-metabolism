function r_AKGD = AKG_Dehydrogenase(M_AKG, M_COASH, M_NAD, M_SCOA, M_NADH, Rate)
    V_AKGD_mf = 8.311*1e4*60.0;
    K_AKGD_AKG = 80.0e-3; K_AKGD_COASH = 55.0e-3; K_AKGD_NAD = 21.0e-3;
    K_AKGD_NADH = 4.5e-3; K_AKGD_SCOA = 6.9e-3;

    V_AKGD_mf = V_AKGD_mf*Rate;
    
    t0 = (V_AKGD_mf*M_AKG*M_COASH*M_NAD) / (K_AKGD_AKG*K_AKGD_COASH*K_AKGD_NAD);
    t1 = 1.0 + (M_AKG / K_AKGD_AKG) + (M_COASH / K_AKGD_COASH) + (M_NAD / K_AKGD_NAD);
    t2 = ((M_AKG*M_COASH) / (K_AKGD_AKG*K_AKGD_COASH)) + ((M_AKG*M_NAD) / (K_AKGD_AKG*K_AKGD_NAD));
    t3 = ((M_COASH*M_NAD) / (K_AKGD_COASH*K_AKGD_NAD)) + ((M_AKG*M_COASH*M_NAD) / (K_AKGD_AKG*K_AKGD_COASH*K_AKGD_NAD));
    t4 = (M_NADH / K_AKGD_NADH) + (M_SCOA / K_AKGD_SCOA);

    r_AKGD = t0 / (t1 + t2 + t3 + t4);
end