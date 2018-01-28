function f = leftMulMQ(Mblk, Qblk, u)
    % suppose that the matrix on the left side has Mblk and Qblk on two
    % diag lines.
%     n = size(Mblk, 1);
%     f = zeros(n^2, 1);
%     
%     f(1:n) = Mblk * u(1:n) + Qblk * u(n+1:2*n);
%     for j = 2:n-1
%         f(1+(j-1)*n:j*n) = Mblk * u(1+(j-1)*n:j*n) ...
%             + Qblk * (u(1+(j-2)*n:(j-1)*n) + u(1+j*n:(j+1)*n));
%     end
%     f(1+(n-1)*n:n*n) = Mblk * u(1+(n-1)*n:n*n) + Qblk * u(1+(n-2)*n:(n-1)*n);
    
    n = size(Mblk, 1);
    U = reshape(u, n, n);
    MU = Mblk * U;
    QU = Qblk * U;
    MU(:, 1:n-1) = MU(:, 1:n-1) + QU(:, 2:n);
    MU(:, 2:n) = MU(:, 2:n) + QU(:, 1:n-1);
    f = reshape(MU, n^2, 1);
end

