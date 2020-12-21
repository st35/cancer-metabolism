function r_GLUD = Glutamate_Dehydrogenase(M_GLU, M_NAD, M_AKG, M_NADH, M_NH3, Rate)
    V_GLUD_mf = 161.9*1e-3*60.0; V_GLUD_mr = V_GLUD_mf*1e-2;

    K_GLUD_GLU = 12.4; K_GLUD_NAD = 60.7e-3;
    K_GLUD_AKG = 2.0; K_GLUD_NH3 = 13.4; K_GLUD_NADH = 40.0e-3;

    V_GLUD_mf = V_GLUD_mf*Rate;
    V_GLUD_mr = V_GLUD_mr*Rate;

    t0 = ((V_GLUD_mf*M_GLU*M_NAD) / (K_GLUD_GLU*K_GLUD_NAD)) - ((V_GLUD_mr*M_AKG*M_NH3*M_NADH) / (K_GLUD_AKG*K_GLUD_NH3*K_GLUD_NADH));
    t1 = 1.0 + (M_GLU / K_GLUD_GLU) + (M_NAD / K_GLUD_NAD) + ((M_GLU*M_NAD) / (K_GLUD_GLU*K_GLUD_NAD));
    t2 = (M_AKG / K_GLUD_AKG) + (M_NH3 / K_GLUD_NH3) + (M_NADH / K_GLUD_NADH);
    t3 = ((M_AKG*M_NH3) / (K_GLUD_AKG*K_GLUD_NH3)) + ((M_AKG*M_NADH) / (K_GLUD_AKG*K_GLUD_NADH)) + ((M_NH3*M_NADH) / (K_GLUD_NH3*K_GLUD_NADH));
    t4 = ((M_AKG*M_NH3*M_NADH) / (K_GLUD_AKG*K_GLUD_NH3*K_GLUD_NADH));

    r_GLUD = t0 / (t1 + t2 + t3 + t4);
end