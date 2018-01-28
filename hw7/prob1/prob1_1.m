clear; clc

N = 64;

errors = [];
times = [];
epses = [1 1e-3 1e-5];
for eps = epses
    aNm1 = obtain_aNm1(N, eps);
    h = 1 / N;

    u_ground = ground_solution(N);
    f = leftMulAh(aNm1, u_ground) / (h^2);
    u0 = ones((N-1)^2, 1);
    
    start = tic;
    [Lblock, Kblock] = LUDecomp_SpeTri(aNm1, N);
    intermediate = solveLEqn_LKb(Lblock, Kblock, f);
    u = solveUEqn_LKb(Lblock, Kblock, intermediate);
    time = toc(start);
    
    error = norm(u-u_ground) * h;
    errors = [errors error];
    times = [times time];
    fprintf('Error is %.4e (eps=%.1e).\n', error, eps);
end

table(epses', errors', times', 'VariableNames', {'eps', 'error', 'time'})

function aNm1 = obtain_aNm1(N, eps)
    % give non divide h^2 sub-block
    aNm1 = diag(ones(N-1, 1)*2*(1+eps)) ...
            - diag(ones(N-2, 1)*eps, 1) ...
            - diag(ones(N-2, 1)*eps, -1);
end

function [Lblock, Kblock] = LUDecomp_SpeTri(aNm1, N)
    % input: assume that the matrix to be decomped is in tridiagonal,
    % with aNm1 on the diagonal line, and -I on both quasi-diagonal line.
    %
    % output: diag line : L1, L2, ... quasi-diag line : K1', K2' ...
    Lblock = zeros(N-1, (N-1)*(N-1)); % Lblock = [L1 L2 ...]
    Kblock = zeros(N-1, (N-2)*(N-1)); % Kblock = [K1 K2 ...]
    
    L = Cholesky(aNm1);
    Lblock(:, 1:(N-1)) = L;
    for i = 1:N-2
        K = -invL(L);
        Kblock(:, 1+(i-1)*(N-1):i*(N-1)) = K;
        L = Cholesky(aNm1 - K'*K);
        Lblock(:, 1+i*(N-1):(i+1)*(N-1)) = L;
    end
    
    Lblock = Lblock * N;
    Kblock = Kblock * N;
end

function u_ground = ground_solution(N)
    mesh = linspace(0, 1, N+1);
    sin_vector = sin(pi*mesh(2:N));
    u_ground = reshape(sin_vector' * sin_vector, (N-1)^2, 1);
end

function y = solveUEqn_LKb(Lblock, Kblock, b)
    n = size(Lblock, 1);
    y = zeros(n^2, 1);
    
    L = Lblock(:, 1+(n-1)*n:n^2);
    y(1+(n-1)*n:n^2) = solveUEqn(L', b(1+(n-1)*n:n^2));
    for j = n-1:-1:1
        L = Lblock(:, 1+(j-1)*n:j*n);
        K = Kblock(:, 1+(j-1)*n:j*n);
        y(1+(j-1)*n:j*n) = solveUEqn(L', b(1+(j-1)*n:j*n) - K*y(1+j*n:(j+1)*n));
    end
end

function x = solveLEqn_LKb(Lblock, Kblock, y)
    n = size(Lblock, 1);
    x = zeros(n^2, 1);
    
    L = Lblock(:, 1:n);
    x(1:n) = solveLEqn(L, y(1:n));
    for j = 2:n
        L = Lblock(:, 1+(j-1)*n:j*n);
        K = Kblock(:, 1+(j-2)*n:(j-1)*n);
        x(1+(j-1)*n:j*n) = solveLEqn(L, y(1+(j-1)*n:j*n) - K'*x(1+(j-2)*n:(j-1)*n));
    end
end