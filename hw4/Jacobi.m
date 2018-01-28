function [x, history] = Jacobi(A, b)
    Dinv = diag(diag(A).^-1);
    LU = diag(diag(A)) - A;
    [x, history] = iterFramework(Dinv * LU, Dinv*b, 0*diag(A));
end

