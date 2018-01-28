function iL = invL(L)
    n = size(L, 1);
    iL = eye(n);
    for j = 1:n-1
        iL(j, :) = iL(j, :) / L(j, j);
        iL(j+1:n, :) = iL(j+1:n, :) - L(j+1:n, j) * iL(j, :);
    end
    iL(n, :) = iL(n, :) / L(n, n);
end
