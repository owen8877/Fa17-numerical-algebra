clear; clc
warning off all

Ns = [32 40 48 56 64 72 80 88 96];
times = [];
eps = 1e-1;
for N = Ns
    Ah = getA(N, eps);
    Mh = getM(N, eps);

    start = tic;
    [V, D] = eigs(Ah, Mh, 1, 'sa');
    time = toc(start);
    times = [times time];
end

p = polyfit(log(Ns), log(times), 1);
testN = [32, 64, 128, 256, 512, 1024];
prediction = exp(polyval(p, log(testN)));
table(Ns', times')
table(testN', prediction')

function A = getA(N, eps)
    n = N - 1;
    h = 1 / N;
    ai = [repmat(1:n, 1, n) repmat(2:n, 1, n) repmat(1:n-1, 1, n) repmat(1:n, 1, n-1) repmat(1:n, 1, n-1)];
    oi = [offset(n, 0:n-1) offset(n-1, 0:n-1) offset(n-1, 0:n-1) offset(n, 1:n-1) offset(n, 0:n-2)]*n;
    ai = ai + oi;
    aj = [repmat(1:n, 1, n) repmat(1:n-1, 1, n) repmat(2:n, 1, n) repmat(1:n, 1, n-1) repmat(1:n, 1, n-1)];
    oj = [offset(n, 0:n-1) offset(n-1, 0:n-1) offset(n-1, 0:n-1) offset(n, 0:n-2) offset(n, 1:n-1)]*n;
    aj = aj + oj;
    av = [ones(1, n^2)*2*(1+eps) -eps*ones(1, 2*n*(n-1)) -ones(1, 2*n*(n-1))] / (h^2);
    A = sparse(ai, aj, av);
end

function M = getM(N, ~)
    n = N - 1;
    ai = [repmat(1:n, 1, n) repmat(2:n, 1, n) repmat(1:n-1, 1, n) repmat(1:n, 1, n-1) repmat(2:n, 1, n-1) repmat(1:n, 1, n-1) repmat(1:n-1, 1, n-1)];
    oi = [offset(n, 0:n-1) offset(n-1, 0:n-1) offset(n-1, 0:n-1) offset(n, 1:n-1) offset(n-1, 1:n-1) offset(n, 0:n-2) offset(n-1, 0:n-2)]*n;
    ai = ai + oi;
    aj = [repmat(1:n, 1, n) repmat(1:n-1, 1, n) repmat(2:n, 1, n) repmat(1:n, 1, n-1) repmat(1:n-1, 1, n-1) repmat(1:n, 1, n-1) repmat(2:n, 1, n-1)];
    oj = [offset(n, 0:n-1) offset(n-1, 0:n-1) offset(n-1, 0:n-1) offset(n, 0:n-2) offset(n-1, 0:n-2) offset(n, 1:n-1) offset(n-1, 1:n-1)]*n;
    aj = aj + oj;
    av = [ones(1, n^2)/2 ones(1, (6*n-2)*(n-1))/12];
    M = sparse(ai, aj, av);
end

function o = offset(m, range)
    o = reshape(repmat(range, m, 1), m*numel(range), 1)';
end