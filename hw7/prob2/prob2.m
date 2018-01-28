clear; clc
warning off all

Ns = [32 64 128 256];
epses = [1e0 1e-1 1e-3 1e-5 1e-7];
ls = [5, 10, 20, 40, 80];
resultCollection = cell(numel(Ns), 1);
for index1 = 1:numel(Ns)
    N = Ns(index1);
    n = N-1;
    h = 1 / N;
    resultCollection{index1} = cell(numel(epses), 1);
    for index2 = numel(epses):-1:1
        eps = epses(index2);
        Ah = getA(N, eps);
        Mh = getM(N, eps);
        resultCollection{index1}{index2} = cell(numel(ls), 1);
        for index3 = numel(ls):-1:1
            l = ls(index3);
            U0 = ones(n^2, 1) / n;
            start = tic;
            [rou, U, history] = IFKSM(l, U0, Ah, Mh, 500 * sqrt(80/l));
            time = toc(start);
            fprintf('N% 5d eps% 4.1e l% 4d - time% 4.1f\n', N, eps, l, time);
            result.rou = rou;
            result.U = U;
            result.time = time;
            result.history = history;
            resultCollection{index1}{index2}{index3} = result;
            % surf(1:n, 1:n, reshape(U, n, n));
        end
    end
end

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