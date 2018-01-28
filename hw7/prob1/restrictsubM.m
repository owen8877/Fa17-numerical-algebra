function sM = restrictsubM(depth)
    % the restriction sub matrix from 2^depth grid to 2^(depth-1) grid, 
    % so the dim should be ((2^(depth-1)-1), (2^(depth)-1)).
    m = 2^(depth-1)-1; % after restriction
    n = 2^(depth)-1; % before restriction
    
    rIndex = repmat(1:m, 1, 3);
    cIndex = [1:2:n-1 2:2:n 3:2:n];
    v1h16 = ones(1, m) / 16;
    v1h8 = ones(1, m) / 8;
    value = [v1h16 v1h8 v1h16];
    
    sM = sparse(rIndex, cIndex, value);
end

