function r_SLC1A5 = Glutamine_Uptake(E_Q, C_Q, Rate)
    V_SLC1A5_mf = 7.67; V_SLC1A5_mr = 0.767;
    K_SLC1A5_Q = 4.0;

    V_SLC1A5_mf = V_SLC1A5_mf*Rate;
    V_SLC1A5_mr = V_SLC1A5_mr*Rate;

    t0 = (V_SLC1A5_mf*E_Q / K_SLC1A5_Q) - (V_SLC1A5_mr*C_Q / K_SLC1A5_Q);
    t1 = 1.0 + (E_Q / K_SLC1A5_Q) + (C_Q / K_SLC1A5_Q);

    r_SLC1A5 = t0 / t1;
end