load('prob1_23_data.mat')

% print table
% for index1 = 1:numel(Ns)
%     N = Ns(index1);
%     fprintf('M=N=%d case:\n', N);
%     results = resultCollection{index1};
%     
%     toDisplay = cell2mat(cellfun(@(c) cell2mat(cellfun(@(d) [d.time numel(d.history) d.error], c, 'UniformOutput', false)'), results, 'UniformOutput', false));
%     fprintf('eps\t\t');
%     fprintf('%.1e\t\t\t', epses);
%     fprintf('\n');
%     
%     for i = 1:numel(query)
%         q = query{i};
%         fprintf('% 9s\t', q);
%         fprintf('% 6.2f % 5d % 3.1e\t', toDisplay(i, :));
%         fprintf('\n');
%     end
%     
%     fprintf('\n');
%     
%     figure
%     title(['N=' num2str(N)])
%     for index2 = 1:numel(epses)
%         s = subplot(2, 3, index2);
%         for i = 1:numel(query)
%             semilogy([1 resultCollection{index1}{index2}{i}.history])
%             xlabel('iteration')
%             ylabel('error')
%             hold(s, 'on');
%         end
%         legend(s, query, 'Location', 'southwest');
%         title(s, ['\epsilon=' num2str(epses(index2))]);
%     end
% end

% fix eps=1e-5, plot error-N
Es = [];
for index1 = 1:numel(Ns)-1
    N = Ns(index1);
    results = resultCollection{index1};
    E = cellfun(@(c) c.error, results{5}, 'UniformOutput', false);
    Es = [Es cell2mat(E)'];
end

figure
for i = 1:numel(query)
    semilogy(Es(i, :))
    hold on
    xlabel('N');
    ylabel('error');
end
xticks(1:5)
xticklabels({'32', '64', '128', '256', '512'})
legend(query)
