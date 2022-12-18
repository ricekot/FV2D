function Fn = Roe_Flux_2(WL, WR, NORM, AREA)
  % check for proper number of arguments
  if nargin ~= 4
      error('Four inputs required.');
  end
  if nargout ~= 1
      error('One output required.');
  end

  % check that first input argument is size 4x1
  if size(WL, 2) ~= 4
      error('Input WL must be a n_edgesx4 matrix.');
  end

  % get dimensions of the input matrix
  nedges = size(WL, 1);

  % create the output matrix
  Fn = zeros(nedges, 4);
  
  % call the computational routine
  for i = 1:nedges
    wl = WL(i, :).';
    wr = WR(i, :).';
    norm = [NORM(i);
            NORM(nedges+i)];
    
    fn = compute(wl, wr, norm);
    Fn(i, :) = AREA(i)*fn.';
  end
end


function fn = compute(u_l, u_r, norm)
    gamma = 1.4;
  
    % velocities
    v_l = u_l(2:3)/u_l(1);
    v_r = u_r(2:3)/u_r(1);
  
    p_l=(gamma-1)*(u_l(4)- 0.5*u_l(1)*(v_l(1)^2 + v_l(2)^2));
    p_r=(gamma-1)*(u_r(4)- 0.5*u_r(1)*(v_r(1)^2 + v_r(2)^2));
  
    h_l = (u_l(4)+p_l)/u_l(1);
    h_r = (u_r(4)+p_r)/u_r(1);
  
    sq_rho = sqrt(u_r(1)/u_l(1));
    rrho = 1/(sq_rho+1);
  
    um = rrho*(v_l+sq_rho*v_r);
  
    hm = rrho*(h_l+sq_rho*h_r);
  
    usq=sum(0.5*um.^2);
  
    am_sq = (gamma-1)*(hm-usq);
    am  = sqrt(am_sq);
    unm = sum(um.*norm);
  
    % Compute Euler flux (first part)
    rhoun_l = sum(u_l(2:3).*norm);
    rhoun_r = sum(u_r(2:3).*norm);
  
    fn(1) = rhoun_l + rhoun_r;
    fn(2) = rhoun_l*v_l(1) + rhoun_r*v_r(1) + (p_l+p_r)*norm(1);
    fn(3) = rhoun_l*v_l(2) + rhoun_r*v_r(2) + (p_l+p_r)*norm(2);
    fn(4) = rhoun_l*h_l + rhoun_r*h_r;
  
    du = u_r-u_l;
  
    lambda0 = abs(unm);
    lambdaP = abs(unm+am);
    lambdaM = abs(unm-am);
  
    % Entropy fix
    eps = 0.5*(abs(rhoun_l/u_l(1)-rhoun_r/u_r(1))+ abs(sqrt(gamma*p_l/u_l(1))-sqrt(gamma*p_r/u_r(1))));
    if lambda0 < 2*eps
        lambda0 = 0.25*lambda0^2/eps + eps;
    end
    if lambdaP < 2*eps
        lambdaP = 0.25*lambdaP^2/eps + eps;
    end
    if lambdaM < 2*eps
        lambdaM = 0.25*lambdaM*lambdaM/eps + eps;
    end
  
    a2 = 0.5*(lambdaP+lambdaM)-lambda0;
    a3 = 0.5*(lambdaP-lambdaM)/am;
    a1 = a2*(gamma-1)/am_sq;
    a4 = a3*(gamma-1);

    a5 = usq*du(1)-um(1)*du(2)-um(2)*du(3)+du(4);
    a6 = unm*du(1)-norm(1)*du(2)-norm(2)*du(3);

    aL1 = a1*a5 - a3*a6;
    bL1 = a4*a5 - a2*a6;

    % Compute Euler flux (second part)
    fn = [
        0.5 * (fn(1) - (lambda0*du(1)+aL1));
        0.5 * (fn(2) - (lambda0*du(2)+aL1*um(1)+bL1*norm(1)));
        0.5 * (fn(3) - (lambda0*du(3)+aL1*um(2)+bL1*norm(2)));
        0.5 * (fn(4) - (lambda0*du(4)+aL1*hm+bL1*unm));
    ];
end