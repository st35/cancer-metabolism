function r_PGK = Phosphoglycerate_Kinase(C_13BPG, C_MgADP, C_3PG, C_MgATP)
    V_PGK_mf = 5.96e4; V_PGK_mr = 2.39e4;
    K_PGK_13BPG = 0.002; K_PGK_MgADP = 0.1; K_PGK_3PG = 0.002; K_PGK_MgATP = 1.0;
    K_PGK_13BPG_i = 1.6; K_PGK_MgADP_i = 0.08; K_PGK_3PG_i = 0.205; K_PGK_MgATP_i = 0.186;

    t0 = ((V_PGK_mf*C_13BPG*C_MgADP) / (K_PGK_MgADP_i*K_PGK_13BPG)) - ((V_PGK_mr*C_3PG*C_MgATP) / (K_PGK_MgATP_i*K_PGK_3PG));
    t1 = 1.0 + (C_13BPG / K_PGK_13BPG_i) + (C_MgADP / K_PGK_MgADP_i);
    t2 = (C_13BPG*C_MgADP) / (K_PGK_MgADP_i*K_PGK_13BPG);
    t3 = (C_3PG / K_PGK_3PG_i) + (C_MgATP / K_PGK_MgATP_i);
    t4 = (C_3PG*C_MgATP) / (K_PGK_MgATP_i*K_PGK_3PG);

    r_PGK = t0 / (t1 + t2 + t3 + t4);
end