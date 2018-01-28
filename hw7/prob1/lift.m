function f = lift(r)
    n = numel(r);
    depth = log2(sqrt(n)+1);
    f = full(liftM(depth)*r);
end