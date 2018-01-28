function [c, s] = givens(a, b)
    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            tao = a / b;
            s = 1 / sqrt(1+tao^2);
            c = s * tao;
        else
            tao = b / a;
            c = 1 / sqrt(1+tao^2);
            s = c * tao;
        end
    end
end