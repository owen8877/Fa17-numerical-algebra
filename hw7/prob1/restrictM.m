function M = restrictM(depth)
    % the restriction matrix from 2^depth grid to 2^(depth-1) grid, so the
    % dim should be ((2^(depth-1)-1)^2, (2^(depth)-1)^2).
%     m = 2^(depth-1)-1; % after restriction
%     n = 2^(depth)-1; % before restriction
%     
%     rIndex = repmat(1:m^2, 1, 9);
%     cIndex = zeros(1, 9*m^2);
%     offset = [0 1 2 n n+1 n+2 2*n 2*n+1 2*n+2];
%     for i = 1:m
%         for j = 1:9
%             cIndex((i-1)*m+1+(j-1)*m^2:i*m+(j-1)*m^2) = 1+2*(i-1)*n+offset(j):2:n+2*(i-1)*n+offset(j)-1;
%         end
%     end
%     
%     v1h16 = ones(1, m^2) / 16;
%     v1h8 = ones(1, m^2) / 8;
%     v1h4 = ones(1, m^2) / 4;
%     value = [v1h16 v1h8 v1h16 v1h8 v1h4 v1h8 v1h16 v1h8 v1h16];
%     
%     M = sparse(rIndex, cIndex, value);
    sM = restrictsubM(depth) * 4;
    M = kron(sM, sM);
end

