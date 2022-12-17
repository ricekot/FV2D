function dUdX = calc_grad_dxinv(U, c2ac, c2nac, dXinv)
    % Check number of input and output arguments
    if nargin ~= 4
        error('calc_grad_dxinv: 4 inputs required.');
    end
    if nargout ~= 1
        error('calc_grad_dxinv: One output required.');
    end

    % Get dimensions of input matrices
    [n_cells_tot, n_fields] = size(U);
    [n_cells, max_nc] = size(c2ac);
    n_dims = size(dXinv, 1);

    % Check number of dimensions
    if n_dims ~= 2
        error('Expected n_dims to be 2; got %i instead.', n_dims);
    end

    % Initialize output matrix
    dUdX = zeros(n_fields, n_dims, n_cells_tot);

    % Loop over cells
    for ic = 1:n_cells
        nc = c2nac(ic);
        % Get the dU matrix
        % TODO: Check +1 -1 stuff below
        dU = zeros(n_fields, max_nc*2);
        for j = 1:nc
            ic2 = c2ac(ic, j);
            dU(:, j) = U(ic2, :)' - U(ic, :)';
        end

        % Get dUdX = dXinv*dU
        for f = 1:n_fields
            for k = 1:n_dims
                for j = 1:nc
                    dUdX(f, k, ic) = dUdX(f, k , ic) ...
                                    + dXinv(k, j, ic) ...
                                    * dU(f, j);
                end
            end
        end
    end
end
