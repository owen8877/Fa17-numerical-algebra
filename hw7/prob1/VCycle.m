function [x, cache] = VCycle(Mblk, Qblk, f, smoothLimit, method, cache)
    n0 = size(Mblk, 1);
    N0 = n0 + 1;
    depth = log2(N0);
    b = f;
    
    xCache = cell(depth, 1);
    bCache = cell(depth, 1);
    MCache = cell(depth, 1);
    QCache = cell(depth, 1);
    
    switch method
        case 'line'
            smoother = @linewiseGS;
        case 'symline'
            smoother = @symlinewiseGS;
        case 'symlinexy'
            smoother = @symlinewisexyGS;
        case 'symliney'
            smoother = @symlinewiseyGS;
        case 'point'
            smoother = @pointwiseGS;
        case 'sympoint'
            smoother = @sympointwiseGS;
    end
    
    if size(cache, 1) == 0
        cache = cell(depth, 1);
    end
    
    for l = depth:-1:2
        bCache{l} = b;
        MCache{l} = Mblk;
        QCache{l} = Qblk;
        c = cache{l};
        [x, c] = smoother(zeros((2^l-1)^2, 1), Mblk, Qblk, b, smoothLimit, c);
        cache{l} = c;
        xCache{l} = x;
        
        r = b - leftMulMQ(Mblk, Qblk, x);
        b = restrict(r);
        
        [Mblk, Qblk] = restrictAh(Mblk, Qblk, l);
    end
    
    % on the coarsest grid, the problem is easy to solve.
    x = b / Mblk;
    
    for l = 2:depth
        x = xCache{l} + lift(x);
        c = cache{l};
        [x, c] = smoother(x, MCache{l}, QCache{l}, bCache{l}, smoothLimit, c);
        cache{l} = c;
    end
end

function [M, Q] = restrictAh(Mt, Qt, depth)
    R = restrictsubM(depth);
    L = liftsubM(depth-1);
    M = R * (6*Mt+8*Qt) * L;
    Q = R * (Mt+4*Qt) * L;
end

function [x, c] = symlinewisexyGS(x0, Mblk, Qblk, f, itrLimit, c)
    if ~isfield(c, 'x')
        c.x = struct();
    end
    if ~isfield(c, 'y')
        c.y = struct();
    end
    [x, cx] = symlinewiseGS(x0, Mblk, Qblk, f, itrLimit, c.x);
    c.x = cx;
    [y, cy] = symlinewiseyGS(x0, Mblk, Qblk, f, itrLimit, c.y);
    c.y = cy;
end