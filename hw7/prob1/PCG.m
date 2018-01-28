function [x, history] = PCG(x0, Mblk, Qblk, b, itrLimit, PCfunc)
    x = x0;
    r0 = b - leftMulMQ(Mblk, Qblk, x);
    r = r0;
    cache = {};
    [z, cache] = PCfunc(r, Mblk, Qblk, cache);
    p = z;
    
    itr = 0;
    history = [];
    while true
        rTz = r' * z;
        Ap = leftMulMQ(Mblk, Qblk, p);
        alpha = rTz / (p'*Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        
        history = [history norm(r)/norm(r0)];
        if (norm(r) < norm(r0) * 1e-6) || itr >= itrLimit
            break
        end
        
        [z, ~] = PCfunc(r, Mblk, Qblk, cache);
        beta = r'*z / rTz;
        p = z + beta * p;
        
        itr = itr + 1;
        if mod(itr, 100) == 0
            fprintf('Iteration %d / %.3e.\n', itr, norm(r)/norm(r0));
        end
    end
%     semilogy(history);
end

