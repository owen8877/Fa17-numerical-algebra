function sM = liftsubM(depth)
    % the lift sub matrix from 2^depth grid to 2^(depth+1) grid, so the
    % dim should be (2^(depth+1)-1, 2^(depth)-1).
    sM = restrictsubM(depth+1)' * 4;
end

