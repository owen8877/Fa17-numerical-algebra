clc; clear

%
N = 64;
ks = [1 3 6 30];

figure
hold on
for j = 1:4
    k = ks(j);
    u = arrayfun(@(x) sin(x*k*pi/N), 1:N-1)';

    % init
    A = diag(repmat(2, N-1, 1)) ...
        + diag(repmat(-1, N-2, 1), 1) ...
        + diag(repmat(-1, N-2, 1), -1);

    D = diag(diag(A));
    L = tril(-A); U = triu(-A);
    M = (D-L) \ U;

    for l = 1:5
        u = (D-L) \ (U*u);
    end
    
    plot(u);
end

legend('1', '3', '6', '30')