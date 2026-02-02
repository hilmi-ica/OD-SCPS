function [P_fail, P_D_post] = BBNFuncs(flow, speed, dH, dP, deta)
    % Load discretization
    if flow >= 700 && flow <= 1100 && speed >= 2877 && speed <= 3356
        L = 1;
    else
        L = 2;
    end
    
    % Degradation discretization
    H_state = 1 + (dH >= 0.10) + (dH >= 0.30);
    P_state = 1 + (dP >= 0.10) + (dP >= 0.30);
    eta_state = 1 + (deta >= 0.05) + (deta >= 0.15);
    
    % CPTs
    P_D_given_L = [
        0.70 0.25 0.05;
        0.30 0.45 0.25
    ];
    
    P_H = [
        0.85 0.10 0.05;
        0.20 0.60 0.20;
        0.05 0.25 0.70
    ];
    
    P_P = [
        0.80 0.15 0.05;
        0.25 0.55 0.20;
        0.10 0.30 0.60
    ];
    
    P_eta = [
        0.90 0.08 0.02;
        0.20 0.60 0.20;
        0.05 0.30 0.65
    ];
    
    P_F_given_D_L = [
        0.005 0.02;
        0.05  0.15;
        0.20  0.45
    ];
    
    % Bayesian update
    for d = 1:3
        likelihood(d) = P_H(d,H_state) * P_P(d,P_state) * ...
            P_eta(d,eta_state);
    end
    
    P_D_post = P_D_given_L(L,:) .* likelihood;
    P_D_post = P_D_post / sum(P_D_post);
    
    % Failure probability
    P_fail = 0;
    for d = 1:3
        P_fail = P_fail + P_D_post(d) * P_F_given_D_L(d,L);
    end
end
