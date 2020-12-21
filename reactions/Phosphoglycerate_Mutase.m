function r_PGM = Phosphoglycerate_Mutase(C_3PG, C_2PG)
    V_PGM_mf = 4.894e5; V_PGM_mr = 4.395e5;
    K_PGM_3PG = 0.168; K_PGM_2PG = 0.0256;

    t0 = (V_PGM_mf*C_3PG / K_PGM_3PG) - (V_PGM_mr*C_2PG / K_PGM_2PG);
    t1 = 1.0 + (C_3PG / K_PGM_3PG) + (C_2PG / K_PGM_2PG);

    r_PGM = t0 / t1;
end