function Fn = viscous_flux(tauL, tauR, phiL, phiR, NORM, AREA, Cp, Pr, R)
    nedges = size(tauL, 1);
    Fn = zeros(nedges, 4);
    
    for ie = 1:nedges
        phil = phiL(ie,:);
        phir = phiR(ie,:);
        taul = tauL(ie,:);
        taur = tauR(ie,:);
        norm = NORM(ie,:);
        area = AREA(ie);
  
        [FL, GL] = compute(taul, phil, Cp, Pr, R);
        [FR, GR] = compute(taur, phir, Cp, Pr, R);
  
        Fn(ie,:) = 0.5*((FL(:)+FR(:))*norm(1) + (GL(:)+GR(:))*norm(2))*area;
    end
end

function [F, G] = compute(grad_u, phi, Cp, Pr, R)
    du_dx = grad_u(1);
    du_dy = grad_u(2);
    dv_dx = grad_u(3);
    dv_dy = grad_u(4);
    dT_dx = grad_u(5);
    dT_dy = grad_u(6);

    u = phi(2);
    v = phi(3);
    T = phi(4) / (R * phi(1));
    
    if isnan(T)
        error('Error! NaN T in viscous_flux.cpp!\n');
    end

    % from HiFiLES: double T_ref = 291.15; double mu_ref = 1.827E-5; 
    % From aerojet.ucdavis.edu/fluenthelp/html/ug/node337.htm:
    T_s = 110.56;      % K - Sutherland's ref. temperature
    T_ref = 273.11;    % K
    mu_ref = 1.716E-5; % kg/m-s

    mu = mu_ref * (T / T_ref)^1.5 * ((T_ref + T_s) / (T + T_s));
    
    diag = (du_dx + dv_dy) / 3.0;

    tauxx = 2.0 * mu * (du_dx - diag);
    tauxy = mu * (du_dy + dv_dx);
    tauyy = 2.0 * mu * (dv_dy - diag);
    Qx = mu * Cp / Pr * dT_dx;
    Qy = mu * Cp / Pr * dT_dy;

    F = [0, -tauxx, -tauxy, -u * tauxx - v * tauxy - Qx];
    G = [0.0, -tauxy, -tauyy, -u * tauxy - v * tauyy - Qy];
end