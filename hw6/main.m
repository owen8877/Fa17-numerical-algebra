clc; clear

H = [2 3 4 5 6; ...
     4 4 5 6 7; ...
     0 3 6 7 8; ...
     0 0 2 8 9; ...
     0 0 0 1 10];
n = 5;
m = n;
itr = 1;

fprintf('Itr        H(2, 1)   H(3, 2)   H(4, 3)   H(5, 4)\n')

while true && itr < 20
    for i = 2:n
        % if abs(H(i, i+1)) < (abs(H(i, i)) + abs(H(i+1, i+1))) * 1e-10
        if abs(H(i, i-1)) < 1e-24
            H(i, i-1) = 0;
        end
    end
    
    % first check for 1x1 block
    while true
        if m >= 2 && H(m, m-1) == 0
            m = m-1;
        % and the 2x2 blocks
        elseif m >= 3 && H(m-1, m-2) == 0
            m = m-2;
        else
            break
        end
    end
    
    if m <= 1
        break
    end
    
    l = m;
    while l >= 2 && H(l, l-1) ~= 0
        l = l-1;
    end
    
    H22 = H(l:m, l:m);
    [Hp, P] = FrancisQR(H22);
    H(1:l-1, l:m) = H(1:l-1, l:m)*P;
    H(l:m, m+1:n) = P'*H(l:m, m+1:n);
    H(l:m, l:m) = Hp;
    
    fprintf('% 2d\t', itr);
    fprintf('% 10.2e', diag(H, -1));
    fprintf('\n');
    itr = itr + 1;
end