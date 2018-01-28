clc
clear

epss = [1, 0.1, 0.01, 0.0001];
a = 0.5;
n = 100;
h = 1/n;

methods = {@SOR, @naiveSOR, @GS, @Jacobi};
names = {'SOR', 'SOR 1.3', 'GS', 'Jacobi'};

figure
for k = 1:4
    eps = epss(k);
    A = zeros(n-1, n-1);
    A(1, 1) = -(2*eps + h);
    for i = 1:n-2
        A(i+1, i+1) = -(2*eps + h);
        A(i, i+1) = eps + h;
        A(i+1, i) = eps;
    end

    b = a * h^2 * ones(n-1, 1);
    b(end) = b(end) - (eps+h);

    yP = zeros(n-1, 1);
    for i = 1:n-1
        x = i/n;
        yP(i) = a * h * i + (1-a) * (1-exp(-x/eps)) / (1-exp(-1/eps));
    end
    
    subplot(2, 2, k)
    fprintf('eps:\t%e\n', eps);
    firstSol = [];
    for j = 1:4
        method = methods{j};
        [x, xhis] = method(A, b);
        if j == 1
            firstSol = x;
        end
        error = historyHandler(xhis, yP);
        semilogy(error)
        ylabel('|y^{(k+1)}-y^{(k)}|')
        xlabel('Iteration')
        hold on
        fprintf('%s\terr:\t%.4e\trel_err:\t%.4e\tstep:\t%d\n', names{j}, norm(yP - x), norm(firstSol - x), size(xhis, 2))
    end
    legend(names)
    title(['\epsilon=' num2str(eps)])
    fprintf('\n')
end

function error = historyHandler(history, yP)
    error = zeros(1, size(history, 2));
    for i = 1:size(history, 2)
        error(i) = norm(history(:, i) - history(:, end));
    end
end