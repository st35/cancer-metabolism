function r_SCOAS = Succinyl_CoA_Synthetase(M_GDP, M_SCOA, M_Pi, M_COASH, M_SUC, M_GTP, Rate)
    V_SCOAS_mf = 1.078*1e4*60.0; V_SCOAS_mr = 2.011*1e4*60.0;
    K_SCOAS_SCOA = 55.0e-3; K_SCOAS_GDP = 8.0e-3; K_SCOAS_Pi = 660.0e-3;
    K_SCOAS_COASH = 20.0e-3; K_SCOAS_GTP = 10.0e-3; K_SCOAS_SUC = 880.0e-3;

    V_SCOAS_mf = V_SCOAS_mf*Rate;
    V_SCOAS_mr = V_SCOAS_mr*Rate;

    t0 = ((V_SCOAS_mf*M_GDP*M_SCOA*M_Pi) / (K_SCOAS_GDP*K_SCOAS_SCOA*K_SCOAS_Pi));
    t1 = ((V_SCOAS_mr*M_COASH*M_SUC*M_GTP) / (K_SCOAS_COASH*K_SCOAS_SUC*K_SCOAS_GTP));

    t2 = 1.0 + (M_GDP / K_SCOAS_GDP) + (M_SCOA / K_SCOAS_SCOA) + (M_Pi / K_SCOAS_Pi);
    t3 = ((M_GDP*M_SCOA) / (K_SCOAS_GDP*K_SCOAS_SCOA));
    t4 = ((M_GDP*M_Pi) / (K_SCOAS_GDP*K_SCOAS_Pi));
    t5 = ((M_SCOA*M_Pi) / (K_SCOAS_SCOA*K_SCOAS_Pi));
    t6 = ((M_GDP*M_SCOA*M_Pi) / (K_SCOAS_GDP*K_SCOAS_SCOA*K_SCOAS_Pi));
    t7 = (M_COASH / K_SCOAS_COASH) + (M_SUC / K_SCOAS_SUC) + (M_GTP / K_SCOAS_GTP);
    t8 = ((M_COASH*M_SUC) / (K_SCOAS_COASH*K_SCOAS_SUC));
    t9 = ((M_COASH*M_GTP) / (K_SCOAS_COASH*K_SCOAS_GTP));
    t10 = ((M_SUC*M_GTP) / (K_SCOAS_SUC*K_SCOAS_GTP));
    t11 = ((M_COASH*M_SUC*M_GTP) / (K_SCOAS_COASH*K_SCOAS_SUC*K_SCOAS_GTP));

    r_SCOAS = (t0 - t1) / (t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 + t11);
end