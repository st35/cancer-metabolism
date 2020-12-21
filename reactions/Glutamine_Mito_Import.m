function r_SLC1A5_var = Glutamine_Mito_Import(C_Q, M_Q, Rate)
    V_SLC1A5_mf = 7.67; V_SLC1A5_mr = 0.767;
    K_SLC1A5_Q = 4.0;

    V_SLC1A5_mf = V_SLC1A5_mf*Rate;
    V_SLC1A5_mr = V_SLC1A5_mr*Rate;

    t0 = (V_SLC1A5_mf*C_Q / K_SLC1A5_Q) - (V_SLC1A5_mr*M_Q / K_SLC1A5_Q);
    t1 = 1.0 + (C_Q / K_SLC1A5_Q) + (M_Q / K_SLC1A5_Q);

    r_SLC1A5_var = t0 / t1;
end