function dydt = model_system(t, y, prm_log)
    
    prm = 10.^prm_log;
    delta       = prm(1, 1);
    base_prd    = prm(2, 1);
    beta        = prm(3, 1);
    k12         = prm(4, 1);
    k23         = prm(5, 1);
    k31         = prm(6, 1);
    
    A = y(1, 1); B = y(2, 1); C = y(3, 1);

    dydt(1, 1)    = - delta * A + base_prd + beta / (1 + (C/k31)^4);
    dydt(2, 1)    = - delta * B + base_prd + beta / (1 + (A/k12)^4);
    dydt(3, 1)    = - delta * C + base_prd + beta / (1 + (B/k23)^4);
    
end