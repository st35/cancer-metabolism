function r_TPI = Triose_Phosphate_Isomerase(C_DHAP, C_GAP)
    V_TPI_mf = 5.10e2; V_TPI_mr = 4.61e1;
    K_TPI_DHAP = 1.62e-1; K_TPI_GAP = 4.30e-1;

    V_TPI_mf = V_TPI_mf*1e2;
    V_TPI_mr = V_TPI_mr*1e2;

    t0 = (V_TPI_mf*C_DHAP / K_TPI_DHAP) - (V_TPI_mr*C_GAP / K_TPI_GAP);
    t1 = 1.0 + (C_DHAP / K_TPI_DHAP) + (C_GAP / K_TPI_GAP);

    r_TPI = t0 / t1;
end