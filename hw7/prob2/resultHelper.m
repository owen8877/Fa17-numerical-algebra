clear; clc

Ns = [32 64 128 256];
epses = [1 1e-1 1e-3 1e-5 1e-7];
ls = [5 10 20 40 80];
lStr = cellfun(@(c) sprintf('l=%d', c), mat2cell(ls, 1, ones(1, numel(ls))), 'UniformOutput', false);

% plot part
figure
for index1 = 1:numel(Ns)
    N = Ns(index1);
    n = N-1;
    h = 1 / N;
    eps = 1e-1;
    l = 80;
    Ah = getA(N, eps);
    Mh = getM(N, eps);
    U0 = ones(n^2, 1) / n;
    start = tic;
    [rou, U, history] = IFKSM(l, U0, Ah, Mh, 500 * sqrt(80/l));
    time = toc(start);
    fprintf('N% 5d eps% 4.1e l% 4d - time% 4.1f\n', N, eps, l, time);
    if U(1, 1) < 0
        U = -U;
    end
    subplot(2, 2, index1);
    mesh((1:n)/N, (1:n)/N, reshape(U, n, n));
    title(sprintf('N=%d \\epsilon=%.1e \\lambda=%.2f', N, eps, rou));
end

% show resultCollection
load resultCollection.mat

for index1 = 1:numel(Ns)
    N = Ns(index1);
    n = N-1;
    fprintf('M=N=%d\n', N);
    fprintf('\t');
    fprintf('eps=%.1e\t\t', epses);
    fprintf('\n');
    for index3 = 1:numel(ls)
        l = ls(index3);
        fprintf('l=% 2d\t', l);
        for index2 = 1:numel(epses)
            eps = epses(index2);
            Mh = getM(N, eps);
            Ugt = reshape(sin((1:n)/n*pi)'*sin((1:n)/n*pi), n^2, 1);
            Ugt = Ugt / sqrt(Ugt'*Mh*Ugt);
            result = resultCollection{index1}{index2}{index3};
            U = result.U;
            if U(1, 1) < 0
                U = -U;
            end
            U = U / sqrt(U'*Mh*U);
            fprintf('% 6.1f % 5d % 4.1e\t', result.time, numel(result.history), norm(Ugt - U)/N);
        end
        fprintf('\n');
    end

    figure
    for index2 = 1:numel(epses)
        eps = epses(index2);
        subplot(2, 3, index2);
        for index3 = 1:numel(ls)
            l = ls(index3);
            result = resultCollection{index1}{index2}{index3};
            semilogy(result.history);
            hold on;
        end
        legend(lStr)
        xlabel('iteration');
        ylabel('relative error');
        title(sprintf('\\epsilon=%.1e', eps))
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