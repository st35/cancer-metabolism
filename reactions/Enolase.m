function r_ENO = Enolase(C_2PG, C_Mg, C_PEP)
    V_ENO_mf = 2.106e4; V_ENO_mr = 5.542e3;
    K_ENO_2PG = 0.046; K_ENO_Mg = 0.14; K_ENO_PEP = 0.11;
    K_ENO_2PG_i = 0.046; K_ENO_Mg_i = 0.14; K_ENO_PEP_i = 0.11;

    t0 = ((V_ENO_mf*C_2PG*C_Mg) / (K_ENO_Mg_i*K_ENO_2PG)) - ((V_ENO_mr*C_PEP*C_Mg) / (K_ENO_Mg_i*K_ENO_PEP));
    t1 = 1.0 + (C_2PG / K_ENO_2PG_i) + (C_Mg / K_ENO_Mg_i) + (C_2PG*C_Mg) / (K_ENO_Mg_i*K_ENO_2PG);
    t2 = (C_PEP / K_ENO_PEP_i) + (C_Mg / K_ENO_Mg_i) + (C_PEP*C_Mg) / (K_ENO_PEP*K_ENO_Mg_i);

    r_ENO = t0 / (t1 + t2);
end