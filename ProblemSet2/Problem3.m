%% Part A
n = 100;
data = 2 * randn(n,1); 

%% Part B
sem = std(data) / sqrt(n);
fprintf('SEM: %.3f\n', sem);

%% Part C
num_bootstrap = 1000;
boot_means = zeros(num_bootstrap,1);

for i = 1:num_bootstrap
    idx = randi(n, n, 1);
    sample = data(idx);
    boot_means(i) = mean(sample);
end

%% Part D
boot_std = std(boot_means);
fprintf('Bootstrap standard deviation: %.3f\n', boot_std);