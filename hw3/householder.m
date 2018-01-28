function [v, b] = householder(x)
    n = numel(x);
    x = x / norm(x, 2);
    sigma = x(2:n)' * x(2:n);
    v = zeros(n, 1);
    v(2:n) = x(2:n);
    if sigma == 0
        b = 0;
    else
        alpha = sqrt(x(1)^2+sigma);
        if x(1) <= 0
            v(1) = x(1) - alpha;
        else
            v(1) = -sigma / (x(1)+alpha);
        end
%         b = 2*v(1)^2/(sigma+v(1)^2);
%         v = v / v(1);
        b = 1;
        v = v / sqrt((sigma+v(1)^2)/2);
    end
end

