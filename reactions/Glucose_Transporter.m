function r_GLUT = Glucose_Transporter(E_GLC, C_GLC, Rate)
    V_GLUT_mf = 7.67; V_GLUT_mr = 0.767;
    K_GLUT_GLC = 5.0e-3;

    V_GLUT_mf = V_GLUT_mf*Rate;
    V_GLUT_mr = V_GLUT_mr*Rate;

    t0 = (V_GLUT_mf*E_GLC) / K_GLUT_GLC;
    t1 = (V_GLUT_mr*C_GLC) / K_GLUT_GLC;
    t2 = 1.0 + (E_GLC / K_GLUT_GLC) + (C_GLC / K_GLUT_GLC);

    r_GLUT = (t0 - t1) / t2;
end