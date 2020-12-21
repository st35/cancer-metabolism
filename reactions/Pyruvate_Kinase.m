function r_PK = Pyruvate_Kinase(C_PEP, C_MgADP, C_PYR, C_MgATP, C_ATP, C_ALA, C_F16BP, C_G16BP)
    V_PK_mf = 2.02e4; V_PK_mr = 47.5;
    K_PK_PEP = 2.25e-1; K_PK_MgADP = 4.74e-1; K_PK_PYR = 4.0; K_PK_MgATP = 3.0;
    K_PK_ATP = 3.39; K_PK_ALA = 0.02;
    K_PK_F16BP = 0.04; K_PK_G16BP = 1.0e-1; 
    L_PK = 0.398;

    t0 = ((V_PK_mf*C_PEP*C_MgADP) / (K_PK_PEP*K_PK_MgADP)) - ((V_PK_mr*C_PYR*C_MgATP) / (K_PK_PYR*K_PK_MgATP));
    t1 = 1.0 + (C_PEP / K_PK_PEP) + (C_MgADP / K_PK_MgADP) + (C_PEP*C_MgADP) / (K_PK_PEP*K_PK_MgADP);
    t2 = (C_PYR / K_PK_PYR) + (C_MgATP / K_PK_MgATP) + (C_PYR*C_MgATP) / (K_PK_PYR*K_PK_MgATP);
    i0 = (1.0 + (C_ATP / K_PK_ATP))^4.0;
    i1 = (1.0 + (C_ALA / K_PK_ALA))^4.0;
    a0 = (1.0 + (C_PEP / K_PK_PEP) + (C_PYR / K_PK_PYR))^4.0;
    a1 = (1.0 + (C_F16BP / K_PK_F16BP) + (C_G16BP / K_PK_G16BP))^4.0;
    N_PK = 1.0 + ((L_PK*i0*i1) / (a0*a1));

    r_PK = (t0) / ((t1 + t2)*N_PK);
end