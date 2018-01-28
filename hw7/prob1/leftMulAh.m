function f = leftMulAh(aNm1, u)
    n = size(aNm1, 1);
    f = leftMulMQ(aNm1, sparse(1:n, 1:n, -ones(1, n)), u);
end

