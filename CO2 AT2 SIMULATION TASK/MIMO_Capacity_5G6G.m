%% =========================================================================
%  MIMO CHANNEL CAPACITY SIMULATION FOR 5G/6G COMMUNICATIONS
%  =========================================================================
%  Simulates capacity under:
%    - Variable antenna configurations (up to Massive MIMO 64x64)
%    - SNR sweep (-10 to 40 dB)
%    - Multiple bandwidths (100 MHz, 400 MHz, 1 GHz)
%    - Fading models: Rayleigh, Rician, CDL-A (5G NR clustered)
%    - Power allocation: Equal & Water-filling
%  =========================================================================

clc; clear; close all;

%% =========================================================================
%  SIMULATION PARAMETERS
%% =========================================================================
params.SNR_dB       = -10:2:40;           % SNR range (dB)
params.SNR_lin      = 10.^(params.SNR_dB/10);
params.N_MC         = 1000;               % Monte Carlo iterations
params.BW_list      = [100e6, 400e6, 1e9]; % Bandwidths: 100MHz, 400MHz, 1GHz
params.BW_labels    = {'100 MHz', '400 MHz', '1 GHz'};
params.ant_configs  = [2, 4, 8, 16, 64]; % Square MIMO: NtxNr
params.K_rician     = [5, 10];           % Rician K-factors
params.fading_types = {'Rayleigh', 'Rician-K5', 'Rician-K10', 'CDL-A'};

fprintf('=============================================================\n');
fprintf('  MIMO Channel Capacity Simulation | 5G/6G Framework\n');
fprintf('=============================================================\n');

%% =========================================================================
%  SECTION 1: CAPACITY vs SNR for Different Antenna Configs (Rayleigh)
%% =========================================================================
fprintf('\n[1/4] Computing Capacity vs SNR for antenna configurations...\n');

figure('Name','Section 1: Capacity vs SNR (Antenna Configs)','Position',[50 50 900 600]);
colors = lines(length(params.ant_configs));
BW = params.BW_list(1); % 100 MHz reference

for a = 1:length(params.ant_configs)
    N = params.ant_configs(a);
    Nt = N; Nr = N;
    C_snr = zeros(1, length(params.SNR_dB));
    
    for s = 1:length(params.SNR_dB)
        rho = params.SNR_lin(s);
        C_mc = zeros(1, params.N_MC);
        
        for mc = 1:params.N_MC
            % Rayleigh fading channel
            H = (randn(Nr,Nt) + 1j*randn(Nr,Nt)) / sqrt(2);
            
            % SVD decomposition
            sigma = svd(H);
            lambda = sigma.^2;  % Eigenvalues
            
            % Water-filling power allocation
            p_wf = water_filling(lambda, rho, Nt);
            
            % Shannon capacity (bps/Hz), scaled by bandwidth
            C_mc(mc) = BW * sum(log2(1 + p_wf .* lambda));
        end
        C_snr(s) = mean(C_mc) / 1e9; % Convert to Gbps
    end
    
    plot(params.SNR_dB, C_snr, '-o', 'Color', colors(a,:), ...
         'LineWidth', 1.8, 'MarkerSize', 4, 'MarkerIndices', 1:3:length(params.SNR_dB));
    hold on;
    
    fprintf('  Antenna %dx%d: Peak = %.2f Gbps @ 40dB SNR\n', N, N, C_snr(end));
end

xlabel('SNR (dB)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Channel Capacity (Gbps)', 'FontSize', 13, 'FontWeight', 'bold');
title({'MIMO Capacity vs SNR — Rayleigh Fading', ...
       sprintf('BW = %s | Water-filling | 5G/6G NR', params.BW_labels{1})}, ...
       'FontSize', 13);
legend(arrayfun(@(x) sprintf('%dx%d MIMO', x, x), params.ant_configs, 'UniformOutput', false), ...
       'Location','northwest', 'FontSize', 10);
grid on; grid minor;
set(gca, 'FontSize', 11);
add_5g_snr_markers(gca);

%% =========================================================================
%  SECTION 2: CAPACITY vs ANTENNA COUNT at Fixed SNRs
%% =========================================================================
fprintf('\n[2/4] Computing Capacity vs Antenna Count...\n');

figure('Name','Section 2: Capacity vs Antenna Count','Position',[100 100 900 600]);
snr_targets = [0, 10, 20, 30];  % Fixed SNR points (dB)
colors2 = [0.2 0.4 0.8; 0.2 0.7 0.3; 0.9 0.5 0.1; 0.8 0.1 0.1];
BW = params.BW_list(2); % 400 MHz

for si = 1:length(snr_targets)
    rho = 10^(snr_targets(si)/10);
    C_ant = zeros(1, length(params.ant_configs));
    
    for a = 1:length(params.ant_configs)
        N = params.ant_configs(a);
        Nt = N; Nr = N;
        C_mc = zeros(1, params.N_MC);
        
        for mc = 1:params.N_MC
            H = (randn(Nr,Nt) + 1j*randn(Nr,Nt)) / sqrt(2);
            sigma = svd(H);
            lambda = sigma.^2;
            p_wf = water_filling(lambda, rho, Nt);
            C_mc(mc) = BW * sum(log2(1 + p_wf .* lambda));
        end
        C_ant(a) = mean(C_mc) / 1e9;
    end
    
    plot(params.ant_configs, C_ant, '-s', 'Color', colors2(si,:), ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', colors2(si,:));
    hold on;
    fprintf('  SNR=%ddB: 2x2=%.1f Gbps → 64x64=%.1f Gbps (%.1fx gain)\n', ...
            snr_targets(si), C_ant(1), C_ant(end), C_ant(end)/max(C_ant(1),0.001));
end

xlabel('Number of Antennas (N_t = N_r)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Channel Capacity (Gbps)', 'FontSize', 13, 'FontWeight', 'bold');
title({'MIMO Capacity vs Antenna Count — Massive MIMO Scaling', ...
       sprintf('BW = %s | Rayleigh | Water-filling', params.BW_labels{2})}, ...
       'FontSize', 13);
legend(arrayfun(@(x) sprintf('SNR = %d dB', x), snr_targets, 'UniformOutput', false), ...
       'Location','northwest', 'FontSize', 11);
set(gca, 'XTick', params.ant_configs);
grid on; grid minor; set(gca, 'FontSize', 11);

% Annotate massive MIMO threshold
xline(16, '--k', '5G Massive MIMO', 'LabelHorizontalAlignment','right', ...
      'FontSize', 10, 'LabelVerticalAlignment','bottom');
xline(64, '--m', '6G Target', 'LabelHorizontalAlignment','left', ...
      'FontSize', 10, 'LabelVerticalAlignment','bottom');

%% =========================================================================
%  SECTION 3: CAPACITY vs FADING MODEL (Rayleigh vs Rician vs CDL-A)
%% =========================================================================
fprintf('\n[3/4] Computing Capacity across Fading Models...\n');

figure('Name','Section 3: Fading Model Comparison','Position',[150 150 900 600]);
Nt = 8; Nr = 8;  % 8x8 MIMO
BW = params.BW_list(2); % 400 MHz
fading_colors = {[0.2 0.4 0.8], [0.1 0.7 0.2], [0.9 0.4 0.1], [0.7 0.1 0.8]};
fading_styles = {'-o', '-s', '-^', '-d'};

for f = 1:length(params.fading_types)
    C_fading = zeros(1, length(params.SNR_dB));
    
    for s = 1:length(params.SNR_dB)
        rho = params.SNR_lin(s);
        C_mc = zeros(1, params.N_MC);
        
        for mc = 1:params.N_MC
            % Generate channel based on fading type
            switch params.fading_types{f}
                case 'Rayleigh'
                    H = generate_rayleigh(Nr, Nt);
                    
                case 'Rician-K5'
                    H = generate_rician(Nr, Nt, 5);
                    
                case 'Rician-K10'
                    H = generate_rician(Nr, Nt, 10);
                    
                case 'CDL-A'
                    H = generate_cdl_a(Nr, Nt);
            end
            
            sigma = svd(H);
            lambda = sigma.^2;
            p_wf = water_filling(lambda, rho, Nt);
            C_mc(mc) = BW * sum(log2(1 + p_wf .* lambda));
        end
        C_fading(s) = mean(C_mc) / 1e9;
    end
    
    plot(params.SNR_dB, C_fading, fading_styles{f}, ...
         'Color', fading_colors{f}, 'LineWidth', 2, ...
         'MarkerSize', 5, 'MarkerIndices', 1:4:length(params.SNR_dB));
    hold on;
    fprintf('  %s: Peak = %.2f Gbps\n', params.fading_types{f}, C_fading(end));
end

xlabel('SNR (dB)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Channel Capacity (Gbps)', 'FontSize', 13, 'FontWeight', 'bold');
title({'MIMO Capacity: Rayleigh vs Rician vs CDL-A Fading', ...
       sprintf('8x8 MIMO | BW = %s | Water-filling', params.BW_labels{2})}, ...
       'FontSize', 13);
legend(params.fading_types, 'Location','northwest', 'FontSize', 11);
grid on; grid minor; set(gca, 'FontSize', 11);
add_5g_snr_markers(gca);

%% =========================================================================
%  SECTION 4: BANDWIDTH & POWER ALLOCATION COMPARISON (Equal vs Water-fill)
%% =========================================================================
fprintf('\n[4/4] Computing Capacity vs Bandwidth & Power Allocation...\n');

figure('Name','Section 4: BW & Power Allocation','Position',[200 200 1100 500]);

% --- Subplot 1: Bandwidth comparison ---
subplot(1,2,1);
Nt = 8; Nr = 8;
bw_colors = {[0.1 0.4 0.9], [0.9 0.4 0.1], [0.1 0.8 0.2]};

for b = 1:length(params.BW_list)
    BW = params.BW_list(b);
    C_bw = zeros(1, length(params.SNR_dB));
    
    for s = 1:length(params.SNR_dB)
        rho = params.SNR_lin(s);
        C_mc = zeros(1, params.N_MC);
        for mc = 1:params.N_MC
            H = generate_rayleigh(Nr, Nt);
            sigma = svd(H);
            lambda = sigma.^2;
            p_wf = water_filling(lambda, rho, Nt);
            C_mc(mc) = BW * sum(log2(1 + p_wf .* lambda));
        end
        C_bw(s) = mean(C_mc) / 1e9;
    end
    
    plot(params.SNR_dB, C_bw, '-', 'Color', bw_colors{b}, 'LineWidth', 2.2);
    hold on;
    fprintf('  BW=%s: Peak = %.2f Gbps\n', params.BW_labels{b}, C_bw(end));
end

xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Capacity (Gbps)', 'FontSize', 12, 'FontWeight', 'bold');
title({'Capacity vs Bandwidth', '8x8 MIMO | Rayleigh | Water-filling'}, 'FontSize', 12);
legend(params.BW_labels, 'Location','northwest', 'FontSize', 10);
grid on; grid minor; set(gca, 'FontSize', 10);

% --- Subplot 2: Equal vs Water-filling ---
subplot(1,2,2);
BW = params.BW_list(2); % 400 MHz
ant_sel = [4, 16]; % Compare 4x4 and 16x16
pwr_colors = {[0.2 0.5 0.9], [0.9 0.2 0.2], [0.1 0.7 0.3], [0.8 0.5 0.0]};
leg_entries = {};
pi_idx = 1;

for a = 1:length(ant_sel)
    N = ant_sel(a);
    Nt = N; Nr = N;
    C_wf = zeros(1, length(params.SNR_dB));
    C_eq = zeros(1, length(params.SNR_dB));
    
    for s = 1:length(params.SNR_dB)
        rho = params.SNR_lin(s);
        C_wf_mc = zeros(1, params.N_MC);
        C_eq_mc = zeros(1, params.N_MC);
        
        for mc = 1:params.N_MC
            H = generate_rayleigh(Nr, Nt);
            sigma = svd(H);
            lambda = sigma.^2;
            rank_H = rank(H);
            
            % Water-filling
            p_wf = water_filling(lambda, rho, Nt);
            C_wf_mc(mc) = BW * sum(log2(1 + p_wf .* lambda));
            
            % Equal power
            p_eq = (rho / rank_H) * ones(rank_H, 1);
            C_eq_mc(mc) = BW * sum(log2(1 + p_eq .* lambda(1:rank_H)));
        end
        C_wf(s) = mean(C_wf_mc) / 1e9;
        C_eq(s) = mean(C_eq_mc) / 1e9;
    end
    
    plot(params.SNR_dB, C_wf, '-', 'Color', pwr_colors{pi_idx}, 'LineWidth', 2.2);
    hold on;
    plot(params.SNR_dB, C_eq, '--', 'Color', pwr_colors{pi_idx+1}, 'LineWidth', 1.8);
    leg_entries{end+1} = sprintf('%dx%d Water-filling', N, N);
    leg_entries{end+1} = sprintf('%dx%d Equal Power', N, N);
    pi_idx = pi_idx + 2;
end

xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Capacity (Gbps)', 'FontSize', 12, 'FontWeight', 'bold');
title({'Water-filling vs Equal Power', sprintf('BW = %s | Rayleigh', params.BW_labels{2})}, 'FontSize', 12);
legend(leg_entries, 'Location','northwest', 'FontSize', 9);
grid on; grid minor; set(gca, 'FontSize', 10);

%% =========================================================================
%  SECTION 5: SUMMARY HEATMAP — Capacity (Gbps) at 20dB SNR
%% =========================================================================
fprintf('\n[Summary] Generating capacity heatmap...\n');

figure('Name','Section 5: Capacity Heatmap','Position',[250 250 800 500]);

BW = params.BW_list(2);
snr_heatmap = [0 5 10 15 20 25 30];
ant_heatmap = [2 4 8 16 32 64];
C_heat = zeros(length(ant_heatmap), length(snr_heatmap));

for a = 1:length(ant_heatmap)
    N = ant_heatmap(a);
    for s = 1:length(snr_heatmap)
        rho = 10^(snr_heatmap(s)/10);
        C_mc = zeros(1, 500);
        for mc = 1:500
            H = generate_rayleigh(N, N);
            sigma = svd(H);
            lambda = sigma.^2;
            p_wf = water_filling(lambda, rho, N);
            C_mc(mc) = BW * sum(log2(1 + p_wf .* lambda));
        end
        C_heat(a,s) = mean(C_mc) / 1e9;
    end
end

imagesc(snr_heatmap, 1:length(ant_heatmap), C_heat);
colormap(turbo);
cb = colorbar;
cb.Label.String = 'Capacity (Gbps)';
cb.Label.FontSize = 12;

for a = 1:length(ant_heatmap)
    for s = 1:length(snr_heatmap)
        text(snr_heatmap(s), a, sprintf('%.1f', C_heat(a,s)), ...
             'HorizontalAlignment','center', 'FontSize', 8.5, ...
             'FontWeight','bold', 'Color','white');
    end
end

yticks(1:length(ant_heatmap));
yticklabels(arrayfun(@(x) sprintf('%dx%d', x, x), ant_heatmap, 'UniformOutput', false));
xlabel('SNR (dB)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Antenna Config (N_t \times N_r)', 'FontSize', 13, 'FontWeight', 'bold');
title({'MIMO Capacity Heatmap (Gbps)', ...
       sprintf('BW = %s | Rayleigh | Water-filling', params.BW_labels{2})}, 'FontSize', 13);
set(gca, 'FontSize', 11);

fprintf('\n=============================================================\n');
fprintf('  Simulation Complete! All plots generated.\n');
fprintf('=============================================================\n');

%% =========================================================================
%  HELPER FUNCTIONS
%% =========================================================================

function H = generate_rayleigh(Nr, Nt)
% Generate normalized Rayleigh MIMO channel
    H = (randn(Nr,Nt) + 1j*randn(Nr,Nt)) / sqrt(2);
end

function H = generate_rician(Nr, Nt, K)
% Generate Rician MIMO channel with K-factor
% LOS component: rank-1 steering vectors
    theta_tx = pi/6;  % AoD
    theta_rx = pi/4;  % AoA
    a_tx = exp(1j*pi*(0:Nt-1)'*sin(theta_tx)) / sqrt(Nt);
    a_rx = exp(1j*pi*(0:Nr-1)'*sin(theta_rx)) / sqrt(Nr);
    H_los = a_rx * a_tx';
    
    % NLOS scatter component
    H_nlos = (randn(Nr,Nt) + 1j*randn(Nr,Nt)) / sqrt(2);
    
    % Combined Rician channel
    H = sqrt(K/(K+1)) * H_los + sqrt(1/(K+1)) * H_nlos;
end

function H = generate_cdl_a(Nr, Nt)
% Simplified 5G NR CDL-A model (Clustered Delay Line)
% Approximates 23 clusters with angular spread
    n_clusters = 23;
    H = zeros(Nr, Nt);
    
    % Cluster powers (dB, normalized from 3GPP CDL-A spec)
    cluster_power_dB = [-13.4 0 -2.2 -4.0 -6.0 -8.2 -9.9 -10.5 ...
                        -7.5 -15.9 -6.6 -16.7 -12.4 -15.2 -10.8 ...
                        -11.3 -12.7 -16.2 -18.3 -18.9 -16.6 -19.9 -29.7];
    cluster_power = 10.^(cluster_power_dB/10);
    cluster_power = cluster_power / sum(cluster_power); % Normalize
    
    % Random AoA/AoD per cluster
    AoA = (2*rand(1,n_clusters)-1) * pi;
    AoD = (2*rand(1,n_clusters)-1) * pi;
    
    for cl = 1:n_clusters
        a_rx = exp(1j*pi*(0:Nr-1)'*sin(AoA(cl))) / sqrt(Nr);
        a_tx = exp(1j*pi*(0:Nt-1)'*sin(AoD(cl))) / sqrt(Nt);
        phase = exp(1j*2*pi*rand());  % Random phase per cluster
        H = H + sqrt(cluster_power(cl)) * phase * (a_rx * a_tx');
    end
    
    % Normalize channel energy
    H = H * sqrt(Nr*Nt) / norm(H,'fro');
end

function p = water_filling(lambda, rho, Nt)
% Water-filling power allocation over MIMO eigenvalues
% lambda: eigenvalues (descending), rho: total SNR, Nt: transmit antennas
    n = length(lambda);
    lambda = sort(lambda, 'descend');
    lambda(lambda < 1e-10) = 1e-10; % Avoid division by zero
    
    % Iterative water-filling
    mu = (rho + sum(1./lambda)) / n;  % Initial water level
    
    for iter = 1:n
        p = max(mu - 1./lambda, 0);
        active = p > 0;
        if sum(active) == 0; break; end
        mu = (rho + sum(1./lambda(active))) / sum(active);
    end
    
    p = max(mu - 1./lambda, 0);
    % Normalize to total power
    if sum(p) > 0
        p = p * (rho / sum(p));
    else
        p = ones(n,1) * rho/n;
    end
end

function add_5g_snr_markers(ax)
% Add 5G/6G reference SNR annotations
    hold(ax, 'on');
    yl = ylim(ax);
    plot(ax, [20 20], yl, ':k', 'LineWidth', 1.2);
    text(ax, 20.5, yl(1) + 0.05*(yl(2)-yl(1)), '5G NR Target', ...
         'FontSize', 8, 'Color', [0.3 0.3 0.3]);
    plot(ax, [35 35], yl, ':m', 'LineWidth', 1.2);
    text(ax, 35.5, yl(1) + 0.05*(yl(2)-yl(1)), '6G Vision', ...
         'FontSize', 8, 'Color', [0.5 0 0.5]);
end
