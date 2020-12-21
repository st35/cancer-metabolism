function r_PFK = Phosphofructokinase(C_F6P, C_MgATP, C_F16BP, C_MgADP, C_ATP, C_Mg, C_23BPG, C_AMP, C_G16BP, C_Pi, C_CIT)
    V_PFK_mf = 15.5e2; V_PFK_mr = 6.78e1;
    K_PFK_F6P = 6.0e-2; K_PFK_MgATP = 6.8e-2; K_PFK_F16BP = 0.65; K_PFK_MgADP = 0.54;
    K_PFK_ATP = 0.1; K_PFK_Mg = 0.2; K_PFK_23BPG = 0.5;
    K_PFK_AMP = 0.3; K_PFK_G16BP = 0.1; K_PFK_Pi = 30.0; K_PFK_CIT = 10.0;
    L_PFK = 2e-3;

    t0 = ((V_PFK_mf*C_F6P*C_MgATP) / (K_PFK_F6P*K_PFK_MgATP)) - ((V_PFK_mr*C_F16BP*C_MgADP) / (K_PFK_F16BP*K_PFK_MgADP));
    t1 = 1.0 + (C_F6P / K_PFK_F6P) + (C_MgATP / K_PFK_MgATP) + ((C_F6P*C_MgATP) / (K_PFK_F6P*K_PFK_MgATP));
    t2 = (C_F16BP / K_PFK_F16BP) + (C_MgADP / K_PFK_MgADP) + ((C_F16BP*C_MgADP) / (K_PFK_F16BP*K_PFK_MgADP));
    i0 = (1.0 + (C_ATP / K_PFK_ATP))^4.0;
    i1 = (1.0 + (C_Mg / K_PFK_Mg))^4.0;
    i2 = (1.0 + (C_23BPG / K_PFK_23BPG))^4.0;
    i3 = (1.0 + (C_CIT / K_PFK_CIT))^4.0;
    a0 = (1.0 + (C_F6P / K_PFK_F6P) + (C_F16BP / K_PFK_F16BP))^4.0;
    a1 = (1.0 + (C_AMP / K_PFK_AMP))^4.0;
    a2 = (1.0 + (C_G16BP / K_PFK_G16BP))^4.0;
    a3 = (1.0 + (C_Pi / K_PFK_Pi))^4.0;
    N_PFK = 1.0 / (1.0 + ((L_PFK*i0*i1*i2*i3) / (a0*a1*a2*a3)));

    r_PFK = N_PFK*((t0) / (t1 + t2));
end