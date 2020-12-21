function r_PIM = Mitochondrial_Pyruvate_Transport(C_PYR, C_H, M_PYR, M_H)
    V_PIM_mf = 6.67e12;

    r_PIM = V_PIM_mf*(C_PYR*C_H - M_PYR*M_H);
end