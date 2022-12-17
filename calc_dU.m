function dU = calc_dU(dUdX, dXE, n_cells)

    % check for proper number of arguments
    if nargin ~= 3
        error('calc_dU: 3 inputs required.');
    end
    if nargout ~= 1
        error('calc_dU: One output required.');
    end
    
    % get dimensions of the input matrices
    [n_fields, n_dims, n_cells_tot] = size(dUdX);
    [~, n_faces, ~] = size(dXE);
    
    % check that first input argument is size 4x1
    if n_dims ~= size(dXE,1)
        error('calc_dU.m: Input matrices must be of compatible sizes!');
    end
    
    if n_cells > n_cells_tot
        error('calc_dU.m: size of dUdX is less than the given n_cells!');
    end
    
    % pre-allocate output matrix
    dU = zeros(n_fields, n_faces, n_cells);

    
    % Calculate delta-U = dU_dX * dX 
    for ic = 1:n_cells
        for face = 1:n_faces
            for field = 1:n_fields
                for dim = 1:n_dims
                    dU(field,face,ic) = ...
                        dU(field,face,ic) ...
                        + dUdX(field,dim,ic) ...
                        * dXE(dim,face,ic);
                end
            end
        end
    end
    
end
    