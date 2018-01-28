function x = multigrid(x0, Mblk, Qblk, f, itrLimit, smoothLimit)
    % suppose that A has Mblk on the main diag line, with Qblk on the
    % quasi-diag line.
    n = size(Mblk, 1);
    r0 = f - leftMulMQ(Mblk, Qblk, x0);
    x = x0;
    
    itr = 0;
    while true
        x = VCycle(Mblk, Qblk, f, smoothLimit, 'symline');
        r = f - leftMulMQ(Mblk, Qblk, x);
        
        itr = itr + 1;
        if mod(itr, 1) == 0
            fprintf('Iteration %d / %.3e.\n', itr, norm(r)/norm(r0));
        end
        
        if (norm(r) < norm(r0) * 6e-5) || itr >= itrLimit
            break
        end
    end
end