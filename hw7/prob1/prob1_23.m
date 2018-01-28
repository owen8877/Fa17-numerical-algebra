clear; clc

query = {'symline' 'symliney' 'sympoint' 'id'};

epses = [1e1 1 1e-1 1e-3 1e-5 1e-7];
Ns = [32 64 128 256 512 1024];
resultCollection = cell(numel(Ns));
% load('prob1_23_data.mat')
for index1 = 1:numel(Ns)
    N = Ns(index1);
    aNm1 = obtain_aNm1(N, eps);
    h = 1 / N;
    aNm1h = aNm1 / (h^2);
    neyeh = -eye(N-1) / (h^2);

    u_ground = ground_solution(N);
    f = leftMulAh(aNm1, u_ground) / (h^2);
    u0 = ones((N-1)^2, 1);
    for index2 = 1:numel(epses)
        eps = epses(index2);
        m = size(query, 2);
        for i = 1:m
            q = query{i};
            fprintf('N %d, eps %.1e, %s\n', N, eps, q);
            
            start = tic;
            if strcmp('id', q)
                [u, history] = PCG(u0, aNm1h, neyeh, f, N, @(r, M, Q, cache) deal(r, cache));
            else
                [u, history] = PCG(u0, sparse(aNm1h), sparse(neyeh), f, 1200/(N/256)^2, ...
                    @(r, M, Q, cache) VCycle(M, Q, r, 1, q, cache));
            end
            result.time = toc(start);
            result.error = norm(u-u_ground) * h;
            result.history = history;
            resultCollection{index1}{index2}{i} = result;
        end
    end
end

save('prob1_23_data.mat', 'Ns', 'epses', 'query', 'resultCollection');

function aNm1 = obtain_aNm1(N, eps)
    % give non divide h^2 sub-block
    aNm1 = diag(ones(N-1, 1)*2*(1+eps)) ...
            - diag(ones(N-2, 1)*eps, 1) ...
            - diag(ones(N-2, 1)*eps, -1);
end

function u_ground = ground_solution(N)
    mesh = linspace(0, 1, N+1);
    sin_vector = sin(pi*mesh(2:N));
    u_ground = reshape(sin_vector' * sin_vector, (N-1)^2, 1);
end