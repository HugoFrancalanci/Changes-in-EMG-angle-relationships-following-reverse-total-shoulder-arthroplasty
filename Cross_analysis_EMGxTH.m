%%  International Shoulder Group Congress 2026 %%

%% Changes in EMG–angle relationships following reverse total shoulder arthroplasty %%
% Hugo Francalanci, Nicolas Holzer, Yosra Cherni, Florent Moissenet

%% ========== PART 1: EMG SCRIPT ========== %%
fprintf('=== PART 1: EMG data generation ===\n');

% addpath(genpath(%%% Define the path to the SPM1D toolbox %%%))

% Define the EMG data files to be analyzed 
data_files = {
    'combined_functional_data_ASYMPTO.mat'
    'combined_functional_data_PREOP.mat'
    'combined_functional_data_POSTOP.mat'
};

% EMG processing function to compute mean activation cycles (tap y on the command window)
plotCombinedEMGPerSubjectWithSPM1DCycles(data_files);

%% ========== PART 2: COMPUTATION OF MEAN EMG CYCLES ========== %%
fprintf('\n=== PART 2: Extraction of mean EMG cycles ===\n');

num_files = length(data_files);
all_data = cell(num_files, 1);
group_labels = {'Asymptomatic', 'Preoperative', 'Postoperative'};

for f = 1:num_files
    data = load(data_files{f});
    all_data{f} = data.individual_data;
end

time_normalized_emg = all_data{1}.time;
muscles_R = all_data{1}.muscles_R;
nb_muscles = length(muscles_R);

muscle_names_clean = cell(nb_muscles, 1);
for m = 1:nb_muscles
    muscle_names_clean{m} = strrep(muscles_R{m}, '_R', '');
    muscle_names_clean{m} = strrep(muscle_names_clean{m}, '_', ' ');
end

% Load previously saved cycle selections to maintain consistency between extracted cycles and results
save_filename = 'cycle_selections';
for f = 1:num_files
    [~, name, ~] = fileparts(data_files{f});
    save_filename = [save_filename '_' name];
end
save_filename = [save_filename '.mat'];

if ~exist(save_filename, 'file')
    error('File of selections not found');
end

loaded_data = load(save_filename);
saved_selections = loaded_data.cycle_selections;
cycle_indices = saved_selections.selections;
normalized_length = 100;
representative_cycles = cell(nb_muscles, num_files);

for m = 1:nb_muscles
    for f = 1:num_files
        data = all_data{f};
        subject_data_R = data.subject_data_R;
        subject_data_L = data.subject_data_L;
        nb_functionals = length(data.functional_labels);
        nb_subjects_R = size(subject_data_R, 1);
        nb_subjects_L = size(subject_data_L, 1);
        all_means_R = zeros(nb_subjects_R, length(time_normalized_emg));
        all_means_L = zeros(nb_subjects_L, length(time_normalized_emg));
        
        for s = 1:nb_subjects_R
            subject_mean = zeros(1, length(time_normalized_emg));
            for func = 1:nb_functionals
                subj_data = subject_data_R{s, func};
                if size(subj_data, 1) == nb_muscles
                    muscle_data = subj_data(m, :);
                else
                    muscle_data = subj_data(:, m)';
                end
                subject_mean = subject_mean + muscle_data;
            end
            all_means_R(s, :) = subject_mean / nb_functionals;
        end
        
        for s = 1:nb_subjects_L
            subject_mean = zeros(1, length(time_normalized_emg));
            for func = 1:nb_functionals
                subj_data = subject_data_L{s, func};
                if size(subj_data, 1) == nb_muscles
                    muscle_data = subj_data(m, :);
                else
                    muscle_data = subj_data(:, m)';
                end
                subject_mean = subject_mean + muscle_data;
            end
            all_means_L(s, :) = subject_mean / nb_functionals;
        end
        
        all_means_combined = [all_means_R; all_means_L];
        saved_indices = cycle_indices{m, f};
        cycles = cell(3, 1);
        
        for c = 1:3
            start_idx = saved_indices(c, 1);
            end_idx = saved_indices(c, 2);
            
            if start_idx > 0 && end_idx > 0
                cycle_range = start_idx:end_idx;
                cycles{c} = all_means_combined(:, cycle_range);
            end
        end
        
        num_subjects = size(all_means_combined, 1);
        single_representative_cycle = zeros(num_subjects, normalized_length);
        
        for s = 1:num_subjects
            normalized_cycles = zeros(3, normalized_length);
            
            for c = 1:3
                if ~isempty(cycles{c}) && size(cycles{c}, 1) >= s
                    cycle_data = cycles{c}(s, :);
                    x_original = linspace(0, 1, length(cycle_data));
                    x_normalized = linspace(0, 1, normalized_length);
                    normalized_cycles(c, :) = interp1(x_original, cycle_data, x_normalized, 'spline');
                end
            end
            
            single_representative_cycle(s, :) = mean(normalized_cycles, 1, 'omitnan');
        end
        
        representative_cycles{m, f} = single_representative_cycle;
    end
end

%% ========== PART 3: LOADING KINEMATIC DATA ========== %%
fprintf('\n=== PART 3: Loading kinematic data ===\n');

% Define the kinematic data files to be analyzed 
load('all_angles_statistics_Asymptomatic.mat');
asymp_kin = all_data;
load('all_angles_statistics_Pre_operatoire.mat');
preop_kin = all_data;
load('all_angles_statistics_Post_operatoire.mat');
postop_kin = all_data;

% Subjects configuration
asymp_subjects = {'A3', 'A4', 'A5', 'A10', 'A12', 'A14', 'A15', 'A17', 'A18', ...
                  'A22', 'A24', 'A29', 'A30', 'A31', 'A32', 'A35', 'A36', 'A37', ...
                  'A38', 'A41'};

% Shoulder side selection
sides_asymp = containers.Map( ...
    {'A3', 'A4', 'A5', 'A10', 'A12', 'A14', 'A15', 'A17', 'A18', ...
     'A22', 'A24', 'A29', 'A30', 'A31', 'A32', 'A35', 'A36', 'A37', ...
     'A38', 'A41'}, ...
    {'left',  'left',  'right',  'right',  'left',  'left',  'left',  'left',  'left', ... 
     'left',  'left',  'left',  'left',  'left',  'left',  'left',  'left',  'left', ...
     'left',  'left'});

preop_subjects = {'S1', 'S3', 'S4', 'S6', 'S7', 'S9', 'S10', 'S11', 'S12', ...
                  'S18', 'S19', 'S20', 'S22', 'S23', 'S26', 'S28', 'S31', 'S32', 'S33', 'S37'};

sides_preop = containers.Map( ...
    {'S1', 'S3', 'S4', 'S6', 'S7', 'S9', 'S10', 'S11', 'S12', ...
     'S18', 'S19', 'S20', 'S22', 'S23', 'S26', 'S28', 'S31', 'S32', 'S33', 'S37'}, ...
    {'right', 'right', 'right', 'right', 'right', 'right', 'left',  'right', 'right', ...
     'left',  'left',  'left',  'left',  'left',  'left',  'left',  'left',  'left', ...
     'right', 'right'});

postop_subjects = preop_subjects;
sides_postop = sides_preop;

% Extraction of thoracohumeral (HT) data by subject
HT_data_asymp = [];
HT_data_preop = [];
HT_data_postop = [];

for i = 1:length(asymp_subjects)
    subj = asymp_subjects{i};
    side = sides_asymp(subj);
    try
        angle = asymp_kin.subjects.(subj).combined.(side).HT.Elevation_globale.mean_deg;
        HT_data_asymp(i, :) = angle(:)';
    catch
        warning('Missing data for %s', subj);
    end
end

for i = 1:length(preop_subjects)
    subj = preop_subjects{i};
    side = sides_preop(subj);
    try
        angle = preop_kin.subjects.(subj).combined.(side).HT.Elevation_globale.mean_deg;
        HT_data_preop(i, :) = angle(:)';
    catch
        warning('Missing HT data for %s', subj);
    end
end

for i = 1:length(postop_subjects)
    subj = postop_subjects{i};
    side = sides_postop(subj);
    try
        angle = postop_kin.subjects.(subj).combined.(side).HT.Elevation_globale.mean_deg;
        HT_data_postop(i, :) = angle(:)';
    catch
        warning('Missing HT data for %s', subj);
    end
end

mean_asymp_kin = mean(HT_data_asymp, 1, 'omitnan');
std_asymp_kin = std(HT_data_asymp, 0, 1, 'omitnan');
mean_preop_kin = mean(HT_data_preop, 1, 'omitnan');
std_preop_kin = std(HT_data_preop, 0, 1, 'omitnan');
mean_postop_kin = mean(HT_data_postop, 1, 'omitnan');
std_postop_kin = std(HT_data_postop, 0, 1, 'omitnan');

%% ========== PART 4: STATISTICAL ANALYSIS ========== %%
fprintf('\n=== PART 4: Individual-level statistical analyses ===\n');

cycle_time_percent = linspace(0, 100, normalized_length);

% Interpolate the kinematics to match the length of EMG cycles
x_kin_original = linspace(0, 100, size(HT_data_asymp, 2));
HT_data_asymp_interp = zeros(size(HT_data_asymp, 1), normalized_length);
HT_data_preop_interp = zeros(size(HT_data_preop, 1), normalized_length);
HT_data_postop_interp = zeros(size(HT_data_postop, 1), normalized_length);

for i = 1:size(HT_data_asymp, 1)
    HT_data_asymp_interp(i, :) = interp1(x_kin_original, HT_data_asymp(i, :), cycle_time_percent, 'spline');
end
for i = 1:size(HT_data_preop, 1)
    HT_data_preop_interp(i, :) = interp1(x_kin_original, HT_data_preop(i, :), cycle_time_percent, 'spline');
end
for i = 1:size(HT_data_postop, 1)
    HT_data_postop_interp(i, :) = interp1(x_kin_original, HT_data_postop(i, :), cycle_time_percent, 'spline');
end

% Data structures for storing results
% Slopes and mean values organized by muscle, group, and subject
HT_by_group = {HT_data_asymp_interp, HT_data_preop_interp, HT_data_postop_interp};
slopes_all = cell(nb_muscles, num_files);
means_all = cell(nb_muscles, num_files);

for m = 1:nb_muscles
    for f = 1:num_files
        emg_data = representative_cycles{m, f};  % (n_sujets x 100)
        ht_data = HT_by_group{f};                % (n_sujets x 100)
        n_subj = size(emg_data, 1);
        slopes = nan(n_subj, 1);
        means = nan(n_subj, 1);
        
        for s = 1:n_subj
            emg_curve = emg_data(s, :);
            ht_curve = ht_data(s, :);
            
            % Linear regression EMG = a * HT + b
            valid_idx = ~isnan(emg_curve) & ~isnan(ht_curve);
            if sum(valid_idx) > 2
                p = polyfit(ht_curve(valid_idx), emg_curve(valid_idx), 1);
                slopes(s) = p(1); % Slope
                means(s) = mean(emg_curve(valid_idx)); % Mean EMG
            end
        end
        
        slopes_all{m, f} = slopes;
        means_all{m, f} = means;
    end
end

%% ========== PART 5: DESCRIPTIVE STATISTICS ========== %%
fprintf('\n=== PART 5: Descriptive statistics ===\n');

% Arrays to store means and standard deviations
slopes_summary = nan(nb_muscles, num_files, 2);  % (:,:,1)=mean, (:,:,2)=SD
means_summary = nan(nb_muscles, num_files, 2);

for m = 1:nb_muscles
    for f = 1:num_files
        slopes_summary(m, f, 1) = mean(slopes_all{m, f}, 'omitnan');
        slopes_summary(m, f, 2) = std(slopes_all{m, f}, 'omitnan');
        means_summary(m, f, 1) = mean(means_all{m, f}, 'omitnan');
        means_summary(m, f, 2) = std(means_all{m, f}, 'omitnan');
    end
end

fprintf('\n--- Slopes ---\n');
fprintf('%-15s', 'Muscle');
for f = 1:num_files
    fprintf('%-25s', group_labels{f});
end
fprintf('\n');

for m = 1:nb_muscles
    fprintf('%-15s', muscles_R{m});
    for f = 1:num_files
        fprintf('%.3f ± %.3f          ', slopes_summary(m, f, 1), slopes_summary(m, f, 2));
    end
    fprintf('\n');
end

fprintf('\n--- Mean EMG ---\n');
fprintf('%-15s', 'Muscle');
for f = 1:num_files
    fprintf('%-25s', group_labels{f});
end
fprintf('\n');

for m = 1:nb_muscles
    fprintf('%-15s', muscles_R{m});
    for f = 1:num_files
        fprintf('%.2f ± %.2f          ', means_summary(m, f, 1), means_summary(m, f, 2));
    end
    fprintf('\n');
end

%% ========== PART 6: ANOVA AND POST-HOC TESTS ========== %%
fprintf('\n=== PART 6: ANOVA and post-hoc tests ===\n');

% Prepare data for ANOVA
% Format: [value, group, subject_id, muscle_id]
slopes_table = [];
means_table = [];

for m = 1:nb_muscles
    for f = 1:num_files
        n = length(slopes_all{m, f});
        slopes_table = [slopes_table; ...
            slopes_all{m, f}, ...
            repmat(f, n, 1), ...
            (1:n)', ...
            repmat(m, n, 1)];  
        means_table = [means_table; ...
            means_all{m, f}, ...
            repmat(f, n, 1), ...
            (1:n)', ...
            repmat(m, n, 1)];
    end
end

slopes_table = slopes_table(~isnan(slopes_table(:,1)), :);
means_table = means_table(~isnan(means_table(:,1)), :);

% ANOVA (slopes)
fprintf('\n--- ANOVA: Slopes ---\n');
[p_slope, tbl_slope, stats_slope] = anovan(slopes_table(:,1), ...
    {slopes_table(:,2), slopes_table(:,4)}, ...
    'model', 'interaction', ...
    'varnames', {'Group', 'Muscle'}, ...
    'display', 'on');

% ANOVA (EMG)
fprintf('\n--- ANOVA: Mean EMG ---\n');
[p_mean, tbl_mean, stats_mean] = anovan(means_table(:,1), ...
    {means_table(:,2), means_table(:,4)}, ...
    'model', 'interaction', ...
    'varnames', {'Group', 'Muscle'}, ...
    'display', 'on');

%% ========== PART 7: POST-HOC TESTS BY MUSCLE ========== %%
fprintf('\n=== PART 7: Post-hoc tests by muscle ===\n');

alpha = 0.05;
posthoc_slopes = cell(nb_muscles, 1);
posthoc_means = cell(nb_muscles, 1);

for m = 1:nb_muscles
    fprintf('\n--- MUSCLE: %s ---\n', muscles_R{m});
    slope_data = [];
    mean_data = [];
    group_id = [];
    
    for f = 1:num_files
        n = length(slopes_all{m, f});
        slope_data = [slope_data; slopes_all{m, f}];
        mean_data = [mean_data; means_all{m, f}];
        group_id = [group_id; repmat(f, n, 1)];
    end
    
    valid = ~isnan(slope_data) & ~isnan(mean_data);
    slope_data = slope_data(valid);
    mean_data = mean_data(valid);
    group_id = group_id(valid);
    
    % Post-hoc tests (slopes)
    fprintf('\nSlopes - Pairwise comparisons (t-tests):\n');
    comparisons_slope = {};
    for i = 1:num_files-1
        for j = i+1:num_files
            data_i = slope_data(group_id == i);
            data_j = slope_data(group_id == j);
            [h, p] = ttest2(data_i, data_j);
            
            sig = '';
            if p < 0.001, sig = '***';
            elseif p < 0.01, sig = '**';
            elseif p < 0.05, sig = '*';
            end
            
            fprintf('  %s vs %s: p=%.4f %s\n', ...
                group_labels{i}, group_labels{j}, p, sig);
            
            comparisons_slope{end+1} = struct('group1', i, 'group2', j, 'p', p, 'sig', h);
        end
    end

    posthoc_slopes{m} = comparisons_slope;
    
    % Post-hoc tests (EMG)
    fprintf('\nEMG - Pairwise comparisons (t-tests):\n');
    comparisons_mean = {};
    for i = 1:num_files-1
        for j = i+1:num_files
            data_i = mean_data(group_id == i);
            data_j = mean_data(group_id == j);
            [h, p] = ttest2(data_i, data_j);
            
            sig = '';
            if p < 0.001, sig = '***';
            elseif p < 0.01, sig = '**';
            elseif p < 0.05, sig = '*';
            end
            
            fprintf('  %s vs %s: p=%.4f %s\n', ...
                group_labels{i}, group_labels{j}, p, sig);
            
            comparisons_mean{end+1} = struct('group1', i, 'group2', j, 'p', p, 'sig', h);
        end
    end
    posthoc_means{m} = comparisons_mean;
end

%% ========== PARTIE 8 : VISUALISATIONS ========== %%
fprintf('\n=== PARTIE 8 : Figure ===\n');

% Figure 1: Comparison of normalized shoulder muscle EMG activity across humeral elevation between groups
colors = {[0, 0, 0], [1, 0, 0], [0, 0, 1]};
line_styles = {'-', '-', '-'};

mean_asymp_kin_interp = mean(HT_data_asymp_interp, 1);
mean_preop_kin_interp = mean(HT_data_preop_interp, 1);
mean_postop_kin_interp = mean(HT_data_postop_interp, 1);

figure('Name', 'Main Figure', 'Color', 'white', 'Position', [50, 50, 1400, 1000]);

for m = 1:nb_muscles
    subplot(3, 2, m);
    hold on;
    box on;
    legend_handles = [];
    legend_labels = {};
    
    for f = 1:num_files
        if ~isempty(representative_cycles{m, f})
            mean_emg = mean(representative_cycles{m, f}, 1, 'omitnan');

            if f == 1
                ht_angle = mean_asymp_kin_interp;
            elseif f == 2
                ht_angle = mean_preop_kin_interp;
            else
                ht_angle = mean_postop_kin_interp;
            end
            
            h_line = plot(ht_angle, mean_emg, line_styles{f}, ...
                'LineWidth', 2, 'Color', colors{f});
            legend_handles = [legend_handles, h_line];
            legend_labels{end+1} = group_labels{f};
        end
    end

    xlabel('Thoracohumeral elevation (°)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Normalised EMG (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title(muscle_names_clean{m}, 'FontSize', 14, 'FontWeight', 'bold');
    xlim([0 80]);
    ylim([0 50]);
    set(gca, 'LineWidth', 1.5, 'FontSize', 11);
    grid off;
    
    if m == 1
        legend(legend_handles, legend_labels, 'Location', 'northeast', 'FontSize', 10);
    end
end

sgtitle('Comparison of normalized shoulder muscle EMG activity across humeral elevation between groups', 'FontSize', 13, 'FontWeight', 'bold');

%% ========== PART 9: RMSD CALCULATION ========== %%
fprintf('\n=== PART 9: RMSD calculation ===\n');

% RMSD between EMG curves of different groups
% For each muscle, compare the mean curves across groups
rmsd_emg_results = zeros(nb_muscles, 3); % Rows: Asymptomatic vs Preoperative, Asymptomatic vs Postoperative, Preoperative vs Postoperative

for m = 1:nb_muscles
    mean_asymp = mean(representative_cycles{m, 1}, 1, 'omitnan');
    mean_preop = mean(representative_cycles{m, 2}, 1, 'omitnan');
    mean_postop = mean(representative_cycles{m, 3}, 1, 'omitnan');
    
    % RMSD (Asymptomatic vs Preoperative)
    rmsd_emg_results(m, 1) = sqrt(mean((mean_asymp - mean_preop).^2, 'omitnan'));
    
    % RMSD (Asymptomatic vs Postoperative)
    rmsd_emg_results(m, 2) = sqrt(mean((mean_asymp - mean_postop).^2, 'omitnan'));
    
    % RMSD (Preoperative vs Postoperative)
    rmsd_emg_results(m, 3) = sqrt(mean((mean_preop - mean_postop).^2, 'omitnan'));
end

% RMSD (EMG)
fprintf('\n--- RMSD between EMG curves ---\n');
fprintf('%-15s %-20s %-20s %-20s\n', 'Muscle', 'Asymp vs Preop', 'Asymp vs Postop', 'Preop vs Postop');
for m = 1:nb_muscles
    fprintf('%-15s %.3f              %.3f              %.3f\n', ...
        muscles_R{m}, rmsd_emg_results(m, 1), rmsd_emg_results(m, 2), rmsd_emg_results(m, 3));
end

% RMSD between kinematic curves (HT)
rmsd_kin = zeros(1, 3);
rmsd_kin(1) = sqrt(mean((mean_asymp_kin_interp - mean_preop_kin_interp).^2, 'omitnan'));
rmsd_kin(2) = sqrt(mean((mean_asymp_kin_interp - mean_postop_kin_interp).^2, 'omitnan'));
rmsd_kin(3) = sqrt(mean((mean_preop_kin_interp - mean_postop_kin_interp).^2, 'omitnan'));

fprintf('\n--- RMSD between kinematic curves (HT) ---\n');
fprintf('Asymp vs Preop:  %.3f°\n', rmsd_kin(1));
fprintf('Asymp vs Postop: %.3f°\n', rmsd_kin(2));
fprintf('Preop vs Postop: %.3f°\n', rmsd_kin(3));

% RMSD by subject (inter-group variability)
rmsd_intra_group = cell(nb_muscles, num_files);

for m = 1:nb_muscles
    for f = 1:num_files
        emg_data = representative_cycles{m, f};
        mean_curve = mean(emg_data, 1, 'omitnan');
        
        n_subj = size(emg_data, 1);
        rmsd_values = zeros(n_subj, 1);
        
        for s = 1:n_subj
            rmsd_values(s) = sqrt(mean((emg_data(s, :) - mean_curve).^2, 'omitnan'));
        end
        
        rmsd_intra_group{m, f} = rmsd_values;
    end
end

% RMSD (inter-group)
fprintf('\n--- RMSD (inter-group) ---\n');
fprintf('%-15s', 'Muscle');
for f = 1:num_files
    fprintf('%-25s', group_labels{f});
end
fprintf('\n');

for m = 1:nb_muscles
    fprintf('%-15s', muscles_R{m});
    for f = 1:num_files
        mean_rmsd = mean(rmsd_intra_group{m, f}, 'omitnan');
        std_rmsd = std(rmsd_intra_group{m, f}, 'omitnan');
        fprintf('%.3f ± %.3f          ', mean_rmsd, std_rmsd);
    end
    fprintf('\n');
end

% RMSD between EMG and kinematic
% Calculate the difference between observed EMG and predicted EMG by linear regression
rmsd_emg_kin_fit = cell(nb_muscles, num_files);

for m = 1:nb_muscles
    for f = 1:num_files
        emg_data = representative_cycles{m, f};
        ht_data = HT_by_group{f};
        n_subj = size(emg_data, 1);
        rmsd_fit = nan(n_subj, 1);
        
        for s = 1:n_subj
            emg_curve = emg_data(s, :);
            ht_curve = ht_data(s, :);
            
            valid_idx = ~isnan(emg_curve) & ~isnan(ht_curve);
            if sum(valid_idx) > 2
                % Linear regression
                p = polyfit(ht_curve(valid_idx), emg_curve(valid_idx), 1);
                emg_predicted = polyval(p, ht_curve(valid_idx));
                rmsd_fit(s) = sqrt(mean((emg_curve(valid_idx) - emg_predicted).^2));
            end
        end
        
        rmsd_emg_kin_fit{m, f} = rmsd_fit;
    end
end

% RMSD (EMG-kinematic fit)
fprintf('\n--- RMSD (EMG-kinematic fit) ---\n');
fprintf('%-15s', 'Muscle');
for f = 1:num_files
    fprintf('%-25s', group_labels{f});
end
fprintf('\n');

for m = 1:nb_muscles
    fprintf('%-15s', muscles_R{m});
    for f = 1:num_files
        mean_rmsd = mean(rmsd_emg_kin_fit{m, f}, 'omitnan');
        std_rmsd = std(rmsd_emg_kin_fit{m, f}, 'omitnan');
        fprintf('%.3f ± %.3f          ', mean_rmsd, std_rmsd);
    end
    fprintf('\n');
end

%% ========== PARTIE 10 : Statistical tests on RMSD ========== %%
fprintf('\n=== PARTIE 10 : Statistical tests on RMSD ===\n');

fprintf('\n--- Inter-group variabilty (ANOVA) ---\n');

alpha = 0.05;
rmsd_intra_table = [];
for m = 1:nb_muscles
    for f = 1:num_files
        n = length(rmsd_intra_group{m, f});
        rmsd_intra_table = [rmsd_intra_table; ...
            rmsd_intra_group{m, f}, ...
            repmat(f, n, 1), ...           % Group
            repmat(m, n, 1)];              % Muscle
    end
end

rmsd_intra_table = rmsd_intra_table(~isnan(rmsd_intra_table(:,1)), :);

[p_intra, tbl_intra, stats_intra] = anovan(rmsd_intra_table(:,1), ...
    {rmsd_intra_table(:,2), rmsd_intra_table(:,3)}, ...
    'model', 'interaction', ...
    'varnames', {'Group', 'Muscle'}, ...
    'display', 'on');

fprintf('\nResults ANOVA inter-group variabilty:\n');
if p_intra(1) < alpha
    fprintf('  Group effect: p = %.4f (significant)\n', p_intra(1));
else
    fprintf('  Group effect: p = %.4f (no significant)\n', p_intra(1));
end

if p_intra(2) < alpha
    fprintf('  Muscle effect: p = %.4f (significant)\n', p_intra(2));
else
    fprintf('  Muscle effect: p = %.4f (no significant)\n', p_intra(2));
end

if p_intra(3) < alpha
    fprintf('  Interaction: p = %.4f (significant)\n', p_intra(3));
else
    fprintf('  Interaction: p = %.4f (no significant)\n', p_intra(3));
end

% Post-hoc tests (if necesssary)
if p_intra(1) < alpha
    fprintf('\n---Post-hoc tests : pairwise comparaison ---\n');
    for i = 1:num_files-1
        for j = i+1:num_files
            data_i = rmsd_intra_table(rmsd_intra_table(:,2) == i, 1);
            data_j = rmsd_intra_table(rmsd_intra_table(:,2) == j, 1);
            
            [h, p] = ttest2(data_i, data_j);
            
            sig = '';
            if p < 0.001
                sig = '***';
            elseif p < 0.01
                sig = '**';
            elseif p < 0.05
                sig = '*';
            end
            
            fprintf('  %s vs %s: p = %.4f %s (Delta = %.3f)\n', ...
                group_labels{i}, group_labels{j}, p, sig, ...
                mean(data_i) - mean(data_j));
        end
    end
end

fprintf('\n--- Fit EMG-HT quality (ANOVA) ---\n');

rmsd_fit_table = [];
for m = 1:nb_muscles
    for f = 1:num_files
        n = length(rmsd_emg_kin_fit{m, f});
        rmsd_fit_table = [rmsd_fit_table; ...
            rmsd_emg_kin_fit{m, f}, ...
            repmat(f, n, 1), ...
            repmat(m, n, 1)];
    end
end

rmsd_fit_table = rmsd_fit_table(~isnan(rmsd_fit_table(:,1)), :);

[p_fit, tbl_fit, stats_fit] = anovan(rmsd_fit_table(:,1), ...
    {rmsd_fit_table(:,2), rmsd_fit_table(:,3)}, ...
    'model', 'interaction', ...
    'varnames', {'Group', 'Muscle'}, ...
    'display', 'on');

fprintf('\nResults ANOVA (fit quality):\n');
if p_fit(1) < alpha
    fprintf('  Group effect: p = %.4f (significant)\n', p_fit(1));
else
    fprintf('  Group effect: p = %.4f (no significant)\n', p_fit(1));
end

if p_fit(2) < alpha
    fprintf('  Muscle effect: p = %.4f (significant)\n', p_fit(2));
else
    fprintf('  Muscle effect: p = %.4f (no significant)\n', p_fit(2));
end

fprintf('\n--- Inter-group variabilty by muscle ---\n');

rmsd_intra_posthoc = cell(nb_muscles, 1);

for m = 1:nb_muscles
    fprintf('\n%s:\n', muscles_R{m});
    data_all = [];
    group_id = [];
    
    for f = 1:num_files
        data_all = [data_all; rmsd_intra_group{m, f}];
        group_id = [group_id; repmat(f, length(rmsd_intra_group{m, f}), 1)];
    end

    valid = ~isnan(data_all);
    data_all = data_all(valid);
    group_id = group_id(valid);
    [p_anova, ~, stats_anova] = anova1(data_all, group_id, 'off');
    
    if p_anova < alpha
        fprintf('  ANOVA: p = %.4f (significant)\n', p_anova);
    else
        fprintf('  ANOVA: p = %.4f (no significant)\n', p_anova);
    end
    
    if p_anova < alpha
        comparisons = {};
        for i = 1:num_files-1
            for j = i+1:num_files
                data_i = data_all(group_id == i);
                data_j = data_all(group_id == j);
                
                [h, p] = ttest2(data_i, data_j);
                
                sig = '';
                if p < 0.001
                    sig = '***';
                elseif p < 0.01
                    sig = '**';
                elseif p < 0.05
                    sig = '*';
                end
                
                fprintf('    %s vs %s: p = %.4f %s\n', ...
                    group_labels{i}, group_labels{j}, p, sig);
                
                comparisons{end+1} = struct('group1', i, 'group2', j, ...
                    'p', p, 'sig', h);
            end
        end
        rmsd_intra_posthoc{m} = comparisons;
    end
end

fprintf('\n--- RMSD EMG between groups (by muscle) ---\n');

rmsd_emg_stats = cell(nb_muscles, 1);

for m = 1:nb_muscles
    fprintf('\n%s:\n', muscles_R{m});
    comparisons = {};
    fprintf('  Asymp vs Preop: RMSD = %.3f%%\n', rmsd_emg_results(m, 1));
    fprintf('  Asymp vs Postop: RMSD = %.3f%%\n', rmsd_emg_results(m, 2));
    fprintf('  Preop vs Postop: RMSD = %.3f%%\n', rmsd_emg_results(m, 3));
    
    for comp = 1:3
        if comp == 1
            data1 = representative_cycles{m, 1};
            data2 = representative_cycles{m, 2};
            label = 'Asymp vs Preop';
        elseif comp == 2
            data1 = representative_cycles{m, 1};
            data3 = representative_cycles{m, 3};
            label = 'Asymp vs Postop';
            data2 = data3;
        else
            data1 = representative_cycles{m, 2};
            data2 = representative_cycles{m, 3};
            label = 'Preop vs Postop';
        end
        
        mean1 = mean(data1, 2, 'omitnan');
        mean2 = mean(data2, 2, 'omitnan');

        if length(mean1) ~= length(mean2)
            [h, p] = ttest2(mean1, mean2);
            test_type = 'independant';
        else
            [h, p] = ttest(mean1, mean2);
            test_type = 'paired';
        end
        
        sig = '';
        if p < 0.001
            sig = '***';
        elseif p < 0.01
            sig = '**';
        elseif p < 0.05
            sig = '*';
        end
        
        fprintf('    %s (test %s): p = %.4f %s\n', label, test_type, p, sig);
    end
end

fprintf('\n--- Correlation RMSD inter-group vs fit quality ---\n');

for f = 1:num_files
    fprintf('\n%s:\n', group_labels{f});
    
    rmsd_intra_vec = [];
    rmsd_fit_vec = [];
    
    for m = 1:nb_muscles
        rmsd_intra_vec = [rmsd_intra_vec; rmsd_intra_group{m, f}];
        rmsd_fit_vec = [rmsd_fit_vec; rmsd_emg_kin_fit{m, f}];
    end
    
    valid = ~isnan(rmsd_intra_vec) & ~isnan(rmsd_fit_vec);
    rmsd_intra_vec = rmsd_intra_vec(valid);
    rmsd_fit_vec = rmsd_fit_vec(valid);
    
    [r, p] = corrcoef(rmsd_intra_vec, rmsd_fit_vec);
    
    if p(1,2) < alpha
        fprintf('  Correlation: r = %.3f, p = %.4f (significant)\n', r(1,2), p(1,2));
    else
        fprintf('  Correlation: r = %.3f, p = %.4f (no significant)\n', r(1,2), p(1,2));
    end
end

