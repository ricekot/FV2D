function Fn = convective_flux(WL, WR, NORM, AREA)
    %CONVECTIVE_FLUX Computes the convective flux at a given interface
    %   between two fluid regions in a 2D CFD simulation
    %
    %   INPUTS:
    %   WL: n_edges x 4 matrix of left state vectors
    %   WR: n_edges x 4 matrix of right state vectors
    %   NORM: n_edges x 2 matrix of interface unit normal vectors
    %   AREA: n_edges x 1 matrix of interface areas
    %
    %   OUTPUTS:
    %   Fn: n_edges x 4 matrix of normal flux vectors

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
        wl = WL(i, :);
        wr = WR(i, :);
        norm = NORM(i, :);
        area = AREA(i);
        Fn(i, :) = compute(wl, wr, norm, area);
    end
end

% u_l and u_r are expected to be row vectors of length 4,
% and norm is expected to be a row vector of length 2.
% The output fn will also be a row vector of length 4.
function fn = compute(u_l, u_r, norm, area)
    gamma = 1.4;

    % velocities
    v_l = u_l(2:3) ./ u_l(1);
    v_r = u_r(2:3) ./ u_r(1);

    % pressure
    p_l = (gamma - 1) * (u_l(4) - 0.5 * u_l(1) * (v_l(1)^2 + v_l(2)^2));
    p_r = (gamma - 1) * (u_r(4) - 0.5 * u_r(1) * (v_r(1)^2 + v_r(2)^2));

    % enthalpy
    h_l = (u_l(4) + p_l) / u_l(1);
    h_r = (u_r(4) + p_r) / u_r(1);

    % face-normal momentum
    rhoun_l = u_l(2:3) * norm';
    rhoun_r = u_r(2:3) * norm';

    % Euler flux
    fn = [
            area * 0.5 * (rhoun_l + rhoun_r);...
            area * 0.5 * (rhoun_l * v_l(1) + rhoun_r * v_r(1) + (p_l + p_r) * norm(1));...
            area * 0.5 * (rhoun_l * v_l(2) + rhoun_r * v_r(2) + (p_l + p_r) * norm(2));...
            area * 0.5 * (rhoun_l * h_l   + rhoun_r * h_r)
        ];
end
