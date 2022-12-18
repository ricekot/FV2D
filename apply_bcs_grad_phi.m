function dphidX = apply_bcs_grad_phi(dphidX, c2c, c2e, bc, unorm, n_cells)
    % Check number of input and output arguments
    if nargin ~= 6
        error('apply_bcs_grad_phi: 6 inputs required.');
    end
    if nargout ~= 1
        error('apply_bcs_grad_phi: One output required.');
    end

    isqrt2 = 1/sqrt(2);

    [n_cells_tot] = size(dphidX);

    for icg = n_cells:n_cells_tot
        ici = c2c(icg);
        ie = c2e(icg); % boundary edge is first edge for ghost cell
        bcg = bc(ie);

        % NOTE: ghost cells ALWAYS on 'right' of all boundary edges
        if bcg == 5
            % Slip Wall Mirror velocities about wall
            % du/dx and dv/dy are the same
            dphidX(1,1,icg) = dphidX(1,1,ici);
            dphidX(2,2,icg) = dphidX(2,2,ici);

            % du/dy and dv/dx are mirrored
            dphidX(1,2,icg) = -dphidX(1,2,ici);
            dphidX(2,1,icg) = -dphidX(2,1,ici);

            % drho/dx and dp/dx are same
            dphidX(1,1,icg) = dphidX(1,1,ici);
            dphidX(4,1,icg) = dphidX(4,1,ici);

            % drho/dy and dp/dy are mirrored
            dphidX(1,2,icg) = -dphidX(1,2,ici);
            dphidX(4,2,icg) = -dphidX(4,2,ici);
        elseif bcg == 6 || bcg == 7
            % Adiabatic Wall or Isothermal Wall
            % Do coordinate transformation at edges from x,y to  normal, tangential coordinates
            n_hat = [unorm(ie,1);  unorm(ie,2)];
            % Find tangential vector
            %t_hat(1) = isqrt2 - (n_hat(1)*isqrt2+n_hat(2)*isqrt2)*n_hat(1)
            %t_hat(2) = isqrt2 - (n_hat(1)*isqrt2+n_hat(2)*isqrt2)*n_hat(2)
            t_hat = cross([isqrt2; isqrt2; 0], [n_hat; 0]);
            t_hat = t_hat(2:3);
            t_hat = t_hat/norm(t_hat);

            % [du/dn; dv/dn]
            dudvdn = dphidX(2:3,:,ici)*n_hat;
            % [du/dt; dv/dt]
            dudvdt = dphidX(2:3,:,ici)*t_hat;
            drhodn = dphidX(1,:,ici)*n_hat;
            dpdn = dphidX(4,:,ici)*n_hat;
            drhodt = dphidX(1,:,ici)*t_hat;
            dpdt = dphidX(4,:,ici)*t_hat;
            
            % dun/dn = dun/du * du/dn + dun/dv * dv/dn, etc.
            % dut/dn, dun/dn are same
            % dut/dt, dun/dt are mirrored
            dun_nt = [dudvdn -dudvdt]'*n_hat;
            dut_nt = [dudvdn -dudvdt]'*t_hat;

            % Transform un, ut to x, y
            dundx = dun_nt'*[n_hat(1); t_hat(1)];
            dundy = dun_nt'*[n_hat(2); t_hat(2)];
            dutdx = dut_nt'*[n_hat(1); t_hat(1)];
            dutdy = dut_nt'*[n_hat(2); t_hat(2)];

            % Transform all derivatives back to x,y
            dphidX(2,1,icg) = dundx/n_hat(1) + dutdx/t_hat(1); % du/dx
            dphidX(2,2,icg) = dundy/n_hat(1) + dutdy/t_hat(1); % du/dy
            dphidX(3,1,icg) = dundx/n_hat(2) + dutdx/t_hat(2); % dv/dx
            dphidX(3,2,icg) = dundy/n_hat(2) + dutdy/t_hat(2); % dv/dy

            % drho/dt, dp/dt are same
            % drho/dn, dp/dn are mirrored
            dphidX(1,1,icg) = -drhodn*n_hat(1) + drhodt*t_hat(1); % drho/dx
            dphidX(1,2,icg) = -drhodn*n_hat(2) + drhodt*t_hat(2); % drho/dy
            dphidX(4,1,icg) = -dpdn*n_hat(1) + dpdt*t_hat(1); % dp/dx
            dphidX(4,2,icg) = -dpdn*n_hat(2) + dpdt*t_hat(2); % dp/dy
        end
    end
end