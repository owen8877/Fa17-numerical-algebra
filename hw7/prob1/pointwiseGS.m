function x = pointwiseGS(x0, Mblk, Qblk, f, itrLimit)
    % suppose that A has Mblk on the main diag line, with Qblk on the
    % quasi-diag line.
    
    % first construct the whole sparse matrix A
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
    
    r0 = f - A*x0;
    x = x0;
    
    g = (D-L) \ f;
    
    itr = 0;
    while true
        x = (D-L) \ (U*x) + g;
        r = f - A*x;
        
        itr = itr + 1;
        if mod(itr, 100) == 0
            fprintf('Iteration %d / %.3e.\n', itr, norm(r)/norm(r0));
        end
        
        if (norm(r) < norm(r0) * 1e-6) || itr >= itrLimit
            break
        end
    end
end