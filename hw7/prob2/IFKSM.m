function [rou, U, history] = IFKSM(l, U0, Ah, Mh, itrLimit)
    rou0 = rayleighQuotient(U0, Ah, Mh);
    rou = rou0;
    r0 = Ah*U0 - rou0*Mh*U0;
    U = U0;
    
    itr = 0;
    history = [];
    while true
        Z = MOBA(U, Ah, Mh, rou, l);
        M = Z' * Mh * Z;
        A = Z' * Ah * Z - rou * M;
        [mu, nu] = minimalEigenPair(A, M);
        rou = rou + mu;
        U = Z * nu;
        
        itr = itr + 1;
        r = Ah*U - rou*Mh*U;
        history = [history norm(r) / norm(r0)];
        if mod(itr, 10) == 0
            fprintf('Itr:% 4d Est error:%.4e\n', itr, norm(r) / norm(r0));
        end
        
        if (norm(r) < norm(r0) * 1e-6) || (itr >= itrLimit)
            break
        end
    end
end