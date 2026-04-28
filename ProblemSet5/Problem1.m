%% Part A
num_trials = length(Kinematics.Trials);
X = cat(1, Kinematics.Trials{:});
[coeff, score, latent, ~, explained] = pca(X);
% coeff (22x22), score (trials concatx22), latent (22x1), explained (22x1)
variance_normalized = explained / 100;
cumulative_variance = cumsum(variance_normalized);
num_pcs = length(variance_normalized);
fprintf('\n--- Variance Explained ---\n');
fprintf('PC1 : %5.2f%%\n', explained(1));
fprintf('PC1 + PC2 : %5.2f%%\n', sum(explained(1:2)));
num_pcs_90 = find(cumulative_variance >= 0.9, 1);
fprintf('PCs for 90%% variance: %d  (cumulative = %.2f%%)\n', num_pcs_90, cumulative_variance(num_pcs_90)*100);


%% ── 5. Figure: scree + cumulative variance ────────────────────────────────
figure('Color','w','Position',[100 100 1100 460]);

% ---- Panel 1: Scree plot (bar) ------------------------------------------
subplot(1,2,1);
barColors = repmat([0.31 0.76 0.97], num_pcs, 1);   % light blue default
barColors(num_pcs_90+1:end,:) = repmat([0.33 0.43 0.48], num_pcs-num_pcs_90, 1);  % grey beyond 90%

b = bar(1:num_pcs, variance_normalized*100, 'FaceColor','flat');
b.CData = barColors;
hold on;

% Annotate first 3 bars
for k = 1:3
    text(k, variance_normalized(k)*100 + 0.4, sprintf('%.1f%%', variance_normalized(k)*100), ...
         'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color',[0.2 0.2 0.2]);
end

xlabel('Principal Component','FontSize',12);
ylabel('Variance Explained (%)','FontSize',12);
title('Scree Plot','FontSize',13,'FontWeight','bold');
xticks(1:num_pcs);
xlim([0.4, num_pcs+0.6]);
box off; grid on; grid minor;

% ---- Panel 2: Cumulative variance ----------------------------------------
subplot(1,2,2);
plot(1:num_pcs, cumulative_variance*100, 'o-', 'Color',[0.31 0.76 0.97], ...
     'LineWidth',2, 'MarkerSize',6, 'MarkerFaceColor','w');
hold on;

% 90% threshold lines
yline(90, '--', '90% threshold', 'Color',[1 0.44 0.26], 'LineWidth',1.5, ...
      'LabelHorizontalAlignment','left','FontSize',9);
xline(num_pcs_90, '--', 'Color',[1 0.44 0.26], 'LineWidth',1.5);
scatter(num_pcs_90, cumulative_variance(num_pcs_90)*100, 80, [1 0.44 0.26], 'filled');
text(num_pcs_90+0.3, cumulative_variance(num_pcs_90)*100 - 3, ...
     sprintf('PC%d: %.1f%%', num_pcs_90, cumulative_variance(num_pcs_90)*100), ...
     'Color',[1 0.44 0.26],'FontSize',9,'FontWeight','bold');

xlabel('Number of PCs','FontSize',12);
ylabel('Cumulative Variance Explained (%)','FontSize',12);
title('Cumulative Variance Explained','FontSize',13,'FontWeight','bold');
xticks(1:num_pcs);
ylim([0 105]);
box off; grid on;

sgtitle('PCA on Kinematics Data — 22 Joint Angles, 337 Trials', ...
        'FontSize',14,'FontWeight','bold');

%% Part B
trial = 1;
cond = 'wr_flexion_l';
cond_ind = find(strcmp(Kinematics.ColumnNames, cond));
temp = Kinematics.Trials{1}(:,cond_ind);