function r_GPI = Glucose_Phosphate_Isomerase(C_G6P, C_F6P)
    V_GPI_mf = 4.8e4; V_GPI_mr = 4.0e4;
    K_GPI_G6P = 0.3; K_GPI_F6P = 0.123;

    K_GPI_G6P = K_GPI_G6P*1e-1;

    t0 = (V_GPI_mf*C_G6P / K_GPI_G6P) - (V_GPI_mr*C_F6P / K_GPI_F6P);
    t1 = 1.0 + (C_G6P / K_GPI_G6P) + (C_F6P / K_GPI_F6P);

    r_GPI = t0 / t1;
end