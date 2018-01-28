function c = powermethod(A)
    criterion = 1e-4;
    n = size(A, 1);
    xo = rand(n, 1); xo = xo / norm(xo);
    xn = A * xo; xn = xn / norm(xn);
    
    y = rand(n, 1); z = rand(n, 1);
    [y, z] = dualNormalize(y, z);
    u = A * y; v = A * z;
    [u, v] = dualNormalize(u, v);
    while true
        xo = xn;
        xn = A * xn;
        xn = xn / norm(xn);
        % prior to check if xn and xo is close enough
        if norm(xn-xo) < criterion
            c = A(1, :)*xn / xn(1);
            return
        end
        
        y = u; z = v;
        zz = z - (y'*z) / (y'*y) * y;
        u = A * y; v = A * z;
        [u, v] = dualNormalize(u, v);
        ur = u - (u'*y)/(y'*y)*y - (u'*zz)/(zz'*zz)*zz;
        vr = v - (v'*y)/(y'*y)*y - (v'*zz)/(zz'*zz)*zz;
        % then check if yz distance close enough
        % disp(norm(ur)+norm(vr))
        if norm(ur)+norm(vr) < criterion
            z = z - (y'*z) / (y'*y) * y;
            y = y / norm(y);
            z = z / norm(z);
            S = [y z]' * A * [y z];
            disp(S)
            c = max(abs(eig(S)));
            return
        end
    end
end

function [yn, zn] = dualNormalize(y, z)
    yn = y / sqrt(y'*y+z'*z);
    zn = z / sqrt(y'*y+z'*z);
end