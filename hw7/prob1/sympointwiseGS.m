function [x, c] = sympointwiseGS(x0, Mblk, Qblk, f, itrLimit, c)
    % suppose that A has Mblk on the main diag line, with Qblk on the
    % quasi-diag line.
    
    % first construct the whole sparse matrix A
    if ~isfield(c, 'A')
        n = size(Mblk, 1);
        [mi, mj, mv] = find(Mblk);
        [qi, qj, qv] = find(Qblk);
        nmblk = numel(mi);
        nqblk = numel(qi);

        av = [repmat(mv', 1, n) repmat(qv', 1, 2*(n-1))];
        ai = [repmat(mi', 1, n) repmat(qi', 1, 2*(n-1))];
        aj = [repmat(mj', 1, n) repmat(qj', 1, 2*(n-1))];
        mshift = reshape(repmat(0:n-1, nmblk, 1), 1, n*nmblk);
        qushift = reshape(repmat(0:n-2, nqblk, 1), 1, (n-1)*nqblk);
        qlshift = reshape(repmat(1:n-1, nqblk, 1), 1, (n-1)*nqblk);
        ai = ai + [mshift qushift qlshift] * n;
        aj = aj + [mshift qlshift qushift] * n;
        A = sparse(ai, aj, av);
        D = diag(diag(A));
        L = -tril(A, -1);
        U = -triu(A, 1);
        
        c.A = A;
        c.D = D;
        c.L = L;
        c.U = U;
    end
    
    r0 = f - c.A*x0;
    x = x0;
    
    g = (c.D-c.L) \ f;
    h = (c.D-c.U) \ f;
    
    itr = 0;
    while true
        x = (c.D-c.L) \ (c.U*x) + g;
        x = (c.D-c.U) \ (c.L*x) + h;
        r = f - c.A*x;
        
        itr = itr + 1;
        if mod(itr, 100) == 0
            fprintf('Iteration %d / %.3e.\n', itr, norm(r)/norm(r0));
        end
        
        if (norm(r) < norm(r0) * 1e-6) || itr >= itrLimit
            break
        end
    end
end