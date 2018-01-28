function M = liftM(depth)
    % the lift matrix from 2^depth grid to 2^(depth+1) grid, so the
    % dim should be ((2^(depth+1)-1)^2, (2^(depth)-1)^2).
    M = restrictM(depth+1)' * 4;
end