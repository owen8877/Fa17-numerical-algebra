function [x, history] = iterFramework(M, g, x0)
    xold = x0;
    xnew = M*x0 + g;
    history = xnew;
    % suppose x0 is quite near to the true solution.
    % q = powermethod(M);
    q = max(abs(eig(M)));
    if q >= 1
        error('Can not handle with matrix norm >= 1.');
    end
    criterion = norm(xnew) / 1e5 * (1-q) / q;
    while true
        xold = xnew;
        xnew = M*xnew + g;
        history = [history xnew];
        if norm(xold-xnew) < criterion
            break
        end
    end
    x = xnew;
end

