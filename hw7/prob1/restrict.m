function f = restrict(r)
    n = numel(r);
    depth = log2(sqrt(n)+1);
    f = full(restrictM(depth)*r);
end