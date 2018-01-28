function Z = MOBA(U, Ah, Mh, rou, l)
    Z = zeros(size(U, 1), l+1);
    Z(:, 1) = U / energyNorm(Mh, U);
    for i = 0:l-1
        w = Ah * Z(:, i+1) - Mh * Z(:, i+1) * rou;
%         ZTMh = Z(:, 1:i+1)' * Mh;
%         for j = 0:i
%             h = Z(:, j+1)' * Mh * w;
%             w = w - h * Z(:, j+1);
%         end
        w = w - Z(:, 1:i+1) * (Z(:, 1:i+1)' * (Mh * w));
        Z(:, i+2) = w / energyNorm(Mh, w);
    end
end

function n = energyNorm(A, x)
    n = sqrt(x'*A*x);
end
