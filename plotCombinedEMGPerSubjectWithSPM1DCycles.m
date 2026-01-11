function plotCombinedEMGPerSubjectWithSPM1DCycles(data_files)
    % 1) This function plots the mean EMG profiles of the four combined
    % movements repeated three times for all subjects, with a separate figure for each muscle,
    % and performs an SPM1D analysis to compare the groups.

    % 2) This function allows selecting three representative cycles,
    % which are then averaged into a SINGLE representative mean cycle

    % Add path to SPM1D toolbox
    addpath(genpath('C:\Users\Francalanci Hugo\Documents\MATLAB\Stage Sainte-Justine\HUG\Statistics\spm1dmatlab-master'))  

    if ischar(data_files) || isstring(data_files)
        data_files = {data_files};
    end

    if length(data_files) > 3
        warning('A maximum of 3 files is supported. Only the first 3 will be processed.');
        data_files = data_files(1:3);
    end

    num_files = length(data_files);
    line_styles = {'-', '-', '-'};
    line_colors = cell(3, 1);
    fill_colors = cell(3, 1);
    default_colors = {
        [0, 0, 0],          % black (Asympto)
        [1, 0, 0],          % red (Preop)
        [0, 0, 1]           % blue (Postop)
    };

    % Colors for confidence intervals (lighter)
    default_fill_colors = {
        [0.5, 0.5, 0.5],    % gray (Asympto)
        [1, 0.5, 0.5],      % light red (Preop)
        [0.5, 0.5, 1]       % light blue (Postop)
    };

    all_data = cell(num_files, 1);
    file_labels = cell(num_files, 1);
    for f = 1:num_files
        if ~exist(data_files{f}, 'file')
            error('The specified data file does not exist: %s', data_files{f});
        end
        fprintf('Loading data from %s...\n', data_files{f});
        data = load(data_files{f});
        if ~isfield(data, 'individual_data')
            error('File %s does not contain the expected "individual_data" structure.', data_files{f});
        end

        all_data{f} = data.individual_data;
        [~, file_name, ~] = fileparts(data_files{f});
        file_labels{f} = file_name;
    end

    group_labels = cell(num_files, 1);
    for f = 1:num_files
        line_colors{f} = default_colors{1};
        fill_colors{f} = default_fill_colors{1};
    end

    if num_files == 1
        group_labels{1} = 'Group';
    elseif num_files == 2
        if contains(lower(file_labels{1}), 'pre')
            group_labels{1} = 'Pre';
            line_colors{1} = default_colors{2};
            fill_colors{1} = default_fill_colors{2};

            if contains(lower(file_labels{2}), 'post')
                group_labels{2} = 'Post';
                line_colors{2} = default_colors{3};
                fill_colors{2} = default_fill_colors{3};
            else
                group_labels{2} = 'Asympt';
                line_colors{2} = default_colors{1}; 
                fill_colors{2} = default_fill_colors{1};
            end
        elseif contains(lower(file_labels{1}), 'asympt')
            group_labels{1} = 'Asympt';
            line_colors{1} = default_colors{1};   % Black for Asympt
            fill_colors{1} = default_fill_colors{1};

            if contains(lower(file_labels{2}), 'pre')
                group_labels{2} = 'Pre';
                line_colors{2} = default_colors{2};  % Red for Pre
                fill_colors{2} = default_fill_colors{2};
            else
                group_labels{2} = 'Post';
                line_colors{2} = default_colors{3};  % Blue for Post
                fill_colors{2} = default_fill_colors{3};
            end
        else
            group_labels{1} = 'Group 1';
            group_labels{2} = 'Group 2';
        end
    elseif num_files == 3
        for f = 1:num_files
            if contains(lower(file_labels{f}), 'asympt')
                group_labels{f} = 'Asympt';
                line_colors{f} = default_colors{1};  
                fill_colors{f} = default_fill_colors{1};
            elseif contains(lower(file_labels{f}), 'pre')
                group_labels{f} = 'Pre';
                line_colors{f} = default_colors{2}; 
                fill_colors{f} = default_fill_colors{2};
            elseif contains(lower(file_labels{f}), 'post')
                group_labels{f} = 'Post';
                line_colors{f} = default_colors{3};  
                fill_colors{f} = default_fill_colors{3};
            else
                group_labels{f} = ['Group ' num2str(f)];
            end
        end
    end

    % Extract necessary information from the first file for initialization
    time_normalized = all_data{1}.time;
    muscles_R = all_data{1}.muscles_R;
    muscles_L = all_data{1}.muscles_L;
    nb_muscles = length(muscles_R);
    if nb_muscles ~= length(muscles_L)
        warning('The number of muscles differs between right and left sides.');
        nb_muscles = min(length(muscles_R), length(muscles_L));
    end

    % Create a structure to store muscle data by file for SPM1D analysis
    all_muscle_data = cell(nb_muscles, 1);
    for m = 1:nb_muscles
        all_muscle_data{m} = cell(num_files, 1);
    end

    % For each muscle, create a separate figure for EMG profiles
    for m = 1:nb_muscles
        figure('Name', sprintf('Muscle %d: %s', m, muscles_R{m}), 'Color', 'white', 'Position', [100+m*50, 100+m*30, 800, 500]);
        hold on;
        legend_handles = [];
        legend_labels = {};
        for f = 1:num_files
            data = all_data{f};
            functional_labels = data.functional_labels;
            subject_data_R = data.subject_data_R;
            subject_data_L = data.subject_data_L;
            subject_ids_R = data.subject_ids_R;
            subject_ids_L = data.subject_ids_L;

            % Number of functional movements
            nb_functionals = length(functional_labels);

            % Number of subjects for each side
            nb_subjects_R = length(subject_ids_R);
            nb_subjects_L = length(subject_ids_L);

            % Matrices to store all subject means
            all_means_R = zeros(nb_subjects_R, length(time_normalized));
            all_means_L = zeros(nb_subjects_L, length(time_normalized));

            % Calculate the mean of 4 movements for each subject on the right side
            for s = 1:nb_subjects_R
                subject_mean = zeros(1, length(time_normalized));

                for func = 1:nb_functionals
                    subj_data = subject_data_R{s, func};

                    if size(subj_data, 1) == nb_muscles
                        % Muscles are in rows
                        muscle_data = subj_data(m, :);
                    else
                        % Muscles are in columns
                        muscle_data = subj_data(:, m)';
                    end

                    subject_mean = subject_mean + muscle_data;
                end

                % Calculate the mean of movements for this subject
                subject_mean = subject_mean / nb_functionals;

                % Store in the array
                all_means_R(s, :) = subject_mean;
            end

            % Calculate the mean of 4 movements for each subject on the left side
            for s = 1:nb_subjects_L
                subject_mean = zeros(1, length(time_normalized));

                for func = 1:nb_functionals
                    subj_data = subject_data_L{s, func};

                    % Check dimensions and extract muscle data
                    if size(subj_data, 1) == nb_muscles
                        % Muscles are in rows
                        muscle_data = subj_data(m, :);
                    else
                        % Muscles are in columns
                        muscle_data = subj_data(:, m)';
                    end

                    subject_mean = subject_mean + muscle_data;
                end

                % Calculate the mean of movements for this subject
                subject_mean = subject_mean / nb_functionals;

                % Store in the array
                all_means_L(s, :) = subject_mean;
            end

            % Combine data from right and left sides
            all_means_combined = [all_means_R; all_means_L];

            % Store data for SPM1D analysis
            all_muscle_data{m}{f} = all_means_combined;

            % Calculate final mean and standard deviation of all subjects
            mean_data = mean(all_means_combined, 1, 'omitnan');
            std_data = std(all_means_combined, 0, 1, 'omitnan');

            % Calculate 95% confidence interval
            n_subjects = size(all_means_combined, 1);
            sem = std_data / sqrt(n_subjects);
            t_crit = tinv(0.975, n_subjects - 1);
            ci_lower = mean_data - t_crit * sem;
            ci_upper = mean_data + t_crit * sem;

            h_line = plot(time_normalized, mean_data, line_styles{f}, 'LineWidth', 2.5, 'Color', line_colors{f});
            x_fill = [time_normalized, fliplr(time_normalized)];
            y_fill = [ci_upper, fliplr(ci_lower)];
            fill(x_fill, y_fill, fill_colors{f}, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            legend_handles = [legend_handles, h_line];

            % Create a label with the file name and number of subjects
            legend_label = sprintf('%s (n=%d)', group_labels{f}, n_subjects);
            legend_labels{end+1} = legend_label;
        end

        % Set the figure title
        if isequal(muscles_R{m}, muscles_L{m})
            title_str = muscles_R{m};
        else
            title_str = [muscles_R{m} ' / ' muscles_L{m}];
        end
        title(title_str, 'FontWeight', 'bold', 'FontSize', 14);

        xlabel('Normalized Time (%)', 'FontSize', 12);
        ylabel('Normalized EMG (% MVC)', 'FontSize', 12);
        xlim([min(time_normalized), max(time_normalized)]);
        grid off;
        box on;

        if ~isempty(legend_handles)
            legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 12);
        end

        set(gca, 'FontSize', 12);
    end

% ========== SAVE/LOAD SELECTION SYSTEM ==========

% Create a save filename based on input files
save_filename = 'cycle_selections';
for f = 1:num_files
    [~, name, ~] = fileparts(data_files{f});
    save_filename = [save_filename '_' name];
end
save_filename = [save_filename '.mat'];

% Variables to store selections
saved_selections = struct();
load_previous = false;

if exist(save_filename, 'file')
    fprintf('\n=== SELECTIONS FILE DETECTED ===\n');
    fprintf('A selections file already exists: %s\n', save_filename);
    
    % Load previous selections
    try
        loaded_data = load(save_filename);
        if isfield(loaded_data, 'cycle_selections')
            saved_selections = loaded_data.cycle_selections;
            
            if isfield(saved_selections, 'nb_muscles') && ...
               isfield(saved_selections, 'num_files') && ...
               isfield(saved_selections, 'muscle_names') && ...
               saved_selections.nb_muscles == nb_muscles && ...
               saved_selections.num_files == num_files
                
                muscles_match = true;
                for m = 1:nb_muscles
                    if ~strcmp(saved_selections.muscle_names{m}, muscles_R{m})
                        muscles_match = false;
                        break;
                    end
                end
                
                if muscles_match
                    fprintf('Compatible selections found\n');
                    fprintf('  - %d muscles\n', saved_selections.nb_muscles);
                    fprintf('  - %d files\n', saved_selections.num_files);
                    fprintf('  - Creation date: %s\n', saved_selections.creation_date);
                    
                    % Ask the user
                    user_choice = input('Do you want to use the previous selections? (y/n): ', 's');
                    if strcmpi(user_choice, 'y') || strcmpi(user_choice, 'yes')
                        load_previous = true;
                        fprintf('→ Loading previous selections...\n');
                    else
                        fprintf('→ New manual selection...\n');
                    end
                else
                    fprintf('Muscle names do not match. New selection required.\n');
                end
            else
                fprintf('Incompatible structure. New selection required.\n');
            end
        end
    catch ME
        fprintf('Error during loading: %s\n', ME.message);
        fprintf('→ New manual selection...\n');
    end
    fprintf('=====================================\n\n');
end

% Structure to store selected cycles and indices
selected_cycles = cell(nb_muscles, num_files);
cycle_indices = cell(nb_muscles, num_files);
representative_cycles = cell(nb_muscles, num_files);
normalized_length = 100;

% Cycle selection loop
for m = 1:nb_muscles
    fprintf('\n===== Cycle selection for muscle %d: %s =====\n', m, muscles_R{m});

    for f = 1:num_files
        % Check whether to load or select
        if load_previous && isfield(saved_selections, 'selections') && ...
           size(saved_selections.selections, 1) >= m && ...
           size(saved_selections.selections, 2) >= f && ...
           ~isempty(saved_selections.selections{m, f})
            
            % Load previous selections
            fprintf('Loading previous selection for %s - %s...\n', muscles_R{m}, group_labels{f});
            
            % Retrieve saved indices
            saved_indices = saved_selections.selections{m, f};
            cycle_indices{m, f} = saved_indices;
            
            % Recreate cycles from indices
            cycles = cell(3, 1);
            for c = 1:3
                start_idx = saved_indices(c, 1);
                end_idx = saved_indices(c, 2);
                
                if start_idx > 0 && end_idx > 0 && start_idx <= end_idx && ...
                   end_idx <= length(time_normalized)
                    
                    cycle_range = start_idx:end_idx;
                    cycles{c} = all_muscle_data{m}{f}(:, cycle_range);
                else
                    fprintf('Invalid indices for cycle %d, using full signal\n', c);
                    cycles{c} = all_muscle_data{m}{f};
                end
            end
            
            selected_cycles{m, f} = cycles;
            
            % Display confirmation graph (optional)
            figure('Name', sprintf('Confirmation - Muscle %d: %s (%s)', m, muscles_R{m}, group_labels{f}), ...
                  'Color', 'white', 'Position', [200+f*50, 200+f*30, 800, 500]);
            
            mean_data = mean(all_muscle_data{m}{f}, 1, 'omitnan');
            plot(time_normalized, mean_data, line_styles{f}, 'LineWidth', 2, 'Color', line_colors{f});
            hold on;
            
            % Display selected cycles
            for c = 1:3
                if ~isempty(cycles{c}) && size(cycles{c}, 2) > 1
                    start_idx = saved_indices(c, 1);
                    end_idx = saved_indices(c, 2);
                    cycle_range = start_idx:end_idx;
                    plot(time_normalized(cycle_range), mean(cycles{c}, 1, 'omitnan'), 'LineWidth', 2, 'Color', [0.2, 0.6, 0.2]);
                end
            end
            
            title(sprintf('Muscle %s - %s: Loaded cycles', muscles_R{m}, group_labels{f}), 'FontSize', 14);
            xlabel('Normalized Time (%)', 'FontSize', 12);
            ylabel('Normalized EMG (% MVC)', 'FontSize', 12);
            legend('Full signal', 'Loaded cycles', 'Location', 'northeast');
            
        else
            % Manual selection
            figure('Name', sprintf('Cycle selection - Muscle %d: %s (%s)', m, muscles_R{m}, group_labels{f}), ...
                  'Color', 'white', 'Position', [200+f*50, 200+f*30, 800, 500]);

            mean_data = mean(all_muscle_data{m}{f}, 1, 'omitnan');
            plot(time_normalized, mean_data, line_styles{f}, 'LineWidth', 2, 'Color', line_colors{f});

            title(sprintf('Muscle %s - %s: Select 3 representative cycles', muscles_R{m}, group_labels{f}), 'FontSize', 14);
            xlabel('Normalized Time (%)', 'FontSize', 12);
            ylabel('Normalized EMG (% MVC)', 'FontSize', 12);
            xlim([min(time_normalized), max(time_normalized)]);
            grid off;

            annotation('textbox', [0.15, 0.01, 0.7, 0.05], 'String', ...
                'Click to select the start and end of each cycle (6 points total)', ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

            fprintf('Please select the start and end points of 3 cycles for %s - %s\n', muscles_R{m}, group_labels{f});
            [x_points, ~] = ginput(6);
            x_points = sort(x_points);

            if length(x_points) ~= 6
                warning('Incorrect number of points selected. Using the full signal.');
                selected_cycles{m, f} = all_muscle_data{m}{f};
                cycle_indices{m, f} = [];  % No valid indices
                continue;
            end

            % Extract cycles and store indices
            cycles = cell(3, 1);
            indices_for_save = zeros(3, 2);  % [start, end] for each cycle
            hold on;

            for c = 1:3
                start_idx = find(time_normalized >= x_points(c*2-1), 1, 'first');
                end_idx = find(time_normalized <= x_points(c*2), 1, 'last');

                if isempty(start_idx) || isempty(end_idx) || start_idx >= end_idx
                    warning('Invalid cycle points. Using the full signal.');
                    selected_cycles{m, f} = all_muscle_data{m}{f};
                    cycle_indices{m, f} = [];
                    continue;
                end

                % Store indices
                indices_for_save(c, 1) = start_idx;
                indices_for_save(c, 2) = end_idx;

                cycle_range = start_idx:end_idx;
                cycles{c} = all_muscle_data{m}{f}(:, cycle_range);

                plot(time_normalized(cycle_range), mean(cycles{c}, 1, 'omitnan'), 'LineWidth', 2, 'Color', [0.2, 0.6, 0.2]);
            end

            selected_cycles{m, f} = cycles;
            cycle_indices{m, f} = indices_for_save;

            title(sprintf('Muscle %s - %s: 3 cycles selected', muscles_R{m}, group_labels{f}), 'FontSize', 14);
            legend('Full signal', 'Selected cycles', 'Location', 'northeast');
        end
        
        % Create representative mean cycle
        num_subjects = size(all_muscle_data{m}{f}, 1);
        single_representative_cycle = zeros(num_subjects, normalized_length);
        
        cycles = selected_cycles{m, f};  % Use cycles (loaded or selected)
        
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

% Save new selections
fprintf('\n=== Saving selections ===\n');
try
    % Prepare save structure
    selections_to_save = struct();
    selections_to_save.nb_muscles = nb_muscles;
    selections_to_save.num_files = num_files;
    selections_to_save.muscle_names = muscles_R;
    selections_to_save.group_labels = group_labels;
    selections_to_save.file_labels = file_labels;
    selections_to_save.creation_date = datestr(now);
    selections_to_save.selections = cycle_indices;  % Store indices
    
    % Save
    cycle_selections = selections_to_save;  % Variable for .mat file
    save(save_filename, 'cycle_selections');
    
    fprintf('Selections saved in: %s\n', save_filename);
    fprintf('  - %d muscles x %d files\n', nb_muscles, num_files);
    fprintf('  - Date: %s\n', datestr(now));
catch ME
    fprintf('Error during saving: %s\n', ME.message);
end
fprintf('=================================\n\n');

    % ========== Figure ==========
    if num_files > 1
        % Create final figure
        figure('Name', 'Article Figure - Representative cycles and SPM1D analyses', ...
               'Color', 'white', 'Position', [50, 50, 1500, 800]);
        
        % Structure to store statistical results
        spm_results = struct();
        
        cycle_time_percent = linspace(0, 100, normalized_length);
        
        for m = 1:nb_muscles

            if isequal(muscles_R{m}, muscles_L{m})
                muscle_name = muscles_R{m};
            else
                muscle_name = [muscles_R{m} '/' muscles_L{m}];
            end
            
            subplot(2, 3, m);
            hold on;
            
            % Plot cycles for each condition
            legend_handles = [];
            legend_labels = {};
            
            for f = 1:num_files
                if ~isempty(representative_cycles{m, f})
                    % Calculate mean and confidence interval
                    mean_data = mean(representative_cycles{m, f}, 1, 'omitnan');
                    std_data = std(representative_cycles{m, f}, 0, 1, 'omitnan');
                    n_subjects = size(representative_cycles{m, f}, 1);
                    
                    sem = std_data / sqrt(n_subjects);
                    t_crit = tinv(0.975, n_subjects - 1);
                    ci_lower = mean_data - t_crit * sem;
                    ci_upper = mean_data + t_crit * sem;
                    h_line = plot(cycle_time_percent, mean_data, line_styles{f}, ...
                        'LineWidth', 2.5, 'Color', line_colors{f});
                    x_fill = [cycle_time_percent, fliplr(cycle_time_percent)];
                    y_fill = [ci_upper, fliplr(ci_lower)];
                    fill(x_fill, y_fill, fill_colors{f}, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
                    legend_handles = [legend_handles, h_line];
                    legend_labels{end+1} = sprintf('%s (n=%d)', group_labels{f}, n_subjects);
                end
            end
            
            % EMG subplot configuration
            title(muscle_name, 'FontWeight', 'bold', 'FontSize', 16);
            xlabel('Cycle (%)', 'FontSize', 14);
            if m == 1 || m == 4
               ylabel('Normalized EMG (%)', 'FontSize', 14);
            else
               set(gca, 'YTickLabel', []); % Remove Y tick labels
               ylabel(''); % Remove Y axis label
            end
            xlim([0, 100]);
            ylim([0, 60]);
            grid off;
            box on;

            if num_files >= 2
            % Define values for each group
            asympt_line = 58.5;
            pre_line = 54.5;
            post_line = 49.5;
    
    has_asympt = any(contains(group_labels, 'Asympt'));
    has_pre = any(contains(group_labels, 'Pre'));
    has_post = any(contains(group_labels, 'Post'));
    
    % Plot lines according to present groups
    if has_asympt
        line([asympt_line, asympt_line], ylim, 'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 2);
    end
    if has_pre
        line([pre_line, pre_line], ylim, 'Color', [1, 0, 0], 'LineStyle', ':', 'LineWidth', 2);
    end
    if has_post
        line([post_line, post_line], ylim, 'Color', [0, 0, 1], 'LineStyle', ':', 'LineWidth', 2);
    end
end
           
% Legend only on first subplot
if m == 1 && ~isempty(legend_handles)
   legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 13);
end

% ========== SPM1D ANALYSIS - BARS UNDER GRAPH ==========

% Initialize results for this muscle
spm_muscle_results = struct();

% Perform pairwise comparisons
comparisons = {};
bar_colors = {};

% ========== ONE-WAY ANOVA - CORRECTED SECTION ==========
% Prepare data for ANOVA
anova_data = {};
anova_labels = {};

for f = 1:num_files
    if ~isempty(representative_cycles{m, f})
        anova_data{end+1} = representative_cycles{m, f};
        anova_labels{end+1} = group_labels{f};
    end
end

% Perform SPM1D ANOVA if at least 2 groups
anova_result = struct();
if length(anova_data) >= 2
    try
        % Check data dimensions
        fprintf('Data dimensions for ANOVA - Muscle %s:\n', muscle_name);
        for i = 1:length(anova_data)
            fprintf('  Group %d (%s): %d subjects x %d time points\n', ...
                i, anova_labels{i}, size(anova_data{i}, 1), size(anova_data{i}, 2));
        end
        
        % Try with data in cell format (standard SPM1D format)
        try
            spm_anova = spm1d.stats.anova1(anova_data);
            fprintf('ANOVA successful with cell format\n');
        catch
            % Convert to matrix format with group vector
            fprintf('Cell format failed, trying with matrix format...\n');
            
            % Combine all data into single matrix
            all_subjects_data = [];
            group_vector = [];
            
            for i = 1:length(anova_data)
                all_subjects_data = [all_subjects_data; anova_data{i}];
                group_vector = [group_vector; repmat(i, size(anova_data{i}, 1), 1)];
            end
            
            % ANOVA call with matrix and group vector
            spm_anova = spm1d.stats.anova1(all_subjects_data, group_vector);
            fprintf('ANOVA successful with matrix format\n');
        end
        
        % Statistical inference
        spmi_anova = spm_anova.inference(0.05, 'interp', true);
        
        % Store ANOVA results
        anova_result.F_stat = spm_anova.z; 
        anova_result.p_values = spmi_anova.p;
        anova_result.significant_clusters = spmi_anova.clusters;
        anova_result.global_significant = ~isempty(spmi_anova.clusters);
        
        fprintf('ANOVA - Muscle %s: ', muscle_name);
        if anova_result.global_significant
        fprintf('Significant effect detected (p < 0.05)\n');
        fprintf('  Number of significant clusters: %d\n', length(spmi_anova.clusters));
        else
        fprintf('No significant effect (p > 0.05)\n');
        end

        catch ME
    fprintf('Error during SPM1D ANOVA for %s: %s\n', muscle_name, ME.message);
    fprintf('Error details: %s\n', ME.getReport);
    anova_result.global_significant = false;
    
    fprintf('Data diagnostic:\n');
    for i = 1:length(anova_data)
        if ~isempty(anova_data{i})
            fprintf('  Group %d: Size = [%d x %d], Min = %.3f, Max = %.3f, NaN = %d\n', ...
                i, size(anova_data{i}, 1), size(anova_data{i}, 2), ...
                min(anova_data{i}(:), [], 'omitnan'), max(anova_data{i}(:), [], 'omitnan'), ...
                sum(isnan(anova_data{i}(:))));
        else
            fprintf('  Group %d: Empty data\n', i);
        end
    end
    end

    else
    fprintf('ANOVA impossible - Less than 2 groups available for %s\n', muscle_name);
    anova_result.global_significant = false;
end

     if num_files == 3
        % All possible comparisons
        comparisons = {'Asympt vs Pre', 'Asympt vs Post', 'Pre vs Post'};
        % Find corresponding indices
        asympt_idx = find(contains(group_labels, 'Asympt'));
        pre_idx = find(contains(group_labels, 'Pre'));
        post_idx = find(contains(group_labels, 'Post'));
        pairs = [asympt_idx pre_idx; asympt_idx post_idx; pre_idx post_idx];
        bar_colors = {[1, 0, 0], [0, 0, 1], [1, 0.5, 0]}; % red, blue, orange
        elseif num_files == 2
               % Single comparison
               comparisons = {sprintf('%s vs %s', group_labels{1}, group_labels{2})};
               pairs = [1 2];
               % Determine color according to comparison type
        if contains(group_labels{1}, 'Asympt') || contains(group_labels{2}, 'Asympt')
               bar_colors = {[0, 0, 0]}; % black
           elseif contains(group_labels{1}, 'Pre') || contains(group_labels{2}, 'Pre')
               bar_colors = {[1, 0, 0]}; % red
           else
               bar_colors = {[0, 0, 1]}; % blue
     end

     end

    % Variables to store p-values
    all_p_values = {};
    % Get current graph limits
    current_ylim = ylim;
    y_min = 0;
    y_max = 60;

% ========== CONDITIONAL T-TESTS ==========
if anova_result.global_significant
   fprintf('Significant ANOVA - Proceeding to post-hoc tests for %s...\n', muscle_name);
   % Perform each comparison
    for comp = 1:size(pairs, 1)
        i = pairs(comp, 1);
        j = pairs(comp, 2);

       if ~isempty(representative_cycles{m, i}) && ~isempty(representative_cycles{m, j})
           data1 = representative_cycles{m, i};
           data2 = representative_cycles{m, j};
    
    % Statistical test
    try
        if size(data1, 1) == size(data2, 1) && ...
           (contains(comparisons{comp}, 'Pre vs Post'))
            % Paired test for Pre vs Post
            spm = spm1d.stats.ttest_paired(data1, data2);
            test_type = 'paired';
        else
            % Independent test
            spm = spm1d.stats.ttest2(data1, data2);
            test_type = 'independent';
        end
        
        % Statistical inference
        spmi = spm.inference(0.05, 'two_tailed', true, 'interp', true);
        
        % Extract significant regions
        if ~isempty(spmi.clusters)
            for cluster = 1:length(spmi.clusters)
                start_pct = (spmi.clusters{cluster}.endpoints(1) - 1) / normalized_length * 100;
                end_pct = (spmi.clusters{cluster}.endpoints(2) - 1) / normalized_length * 100;
                
                % Position of bars under main graph
                bar_height = 1;
                y_position = 5 - (comp * 1); 
                rectangle('Position', [start_pct, y_position, end_pct - start_pct, bar_height], ...
                          'FaceColor', bar_colors{comp}, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            end
            
            % Store p-values (same code as before)
            cluster_info = struct();
            cluster_info.comparison = comparisons{comp};
            cluster_info.test_type = test_type;
            cluster_info.clusters = spmi.clusters;
            cluster_info.p_values = spmi.p;
            
            all_p_values{end+1} = cluster_info;
        end
        
    catch ME
        warning('Error during SPM1D analysis for %s - %s: %s', ...
            muscle_name, comparisons{comp}, ME.message);
    end

    end

    end

       else
           fprintf('Non-significant ANOVA - Post-hoc tests not necessary for %s\n', muscle_name);
           all_p_values = {};
    
end

% store results for this muscle
spm_muscle_results.muscle_name = muscle_name;
spm_muscle_results.anova_result = anova_result;
spm_muscle_results.comparisons = all_p_values;
spm_results.(sprintf('muscle_%d', m)) = spm_muscle_results;

% Add legend for bars on last subplot
    if m == nb_muscles && ~isempty(comparisons)
       legend_handles_spm = [];
        for comp = 1:length(comparisons)
            h_dummy = plot(NaN, NaN, 's', 'MarkerFaceColor', bar_colors{comp}, ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 8);
            legend_handles_spm = [legend_handles_spm, h_dummy];
        end

        legend_combined = [legend_handles, legend_handles_spm];
        legend_labels_combined = [legend_labels, comparisons];
        legend(legend_combined, legend_labels_combined, 'Location', 'best', 'FontSize', 13);
    end

    sgtitle('Representative EMG Cycles and SPM1D Statistical Analyses', 'FontSize', 16, 'FontWeight', 'bold');
        
    % ========== P-VALUES TABLE ==========
    fprintf('\n========== SPM1D ANALYSIS RESULTS ==========\n');
    fprintf('Table of p-values and significant regions:\n\n');
    
            % Table header
            fprintf('%-15s %-20s %-15s %-20s %-15s %-15s\n', ...
                'Muscle', 'Comparison', 'Test Type', 'Cycle Start (%)', 'Cycle End (%)', 'p-value');
            fprintf(repmat('-', 1, 120));
            fprintf('\n');
            
            % Display results for each muscle
            for m = 1:nb_muscles
                muscle_field = sprintf('muscle_%d', m);
                if isfield(spm_results, muscle_field)
                    muscle_results = spm_results.(muscle_field);
                    
                    if ~isempty(muscle_results.comparisons)
                        fprintf('%-15s\n', muscle_results.muscle_name);
                        
                        for comp_idx = 1:length(muscle_results.comparisons)
                            comp_data = muscle_results.comparisons{comp_idx};
                            
                            if ~isempty(comp_data.clusters)
                                for cluster = 1:length(comp_data.clusters)
                                    start_pct = (comp_data.clusters{cluster}.endpoints(1) - 1) / normalized_length * 100;
                                    end_pct = (comp_data.clusters{cluster}.endpoints(2) - 1) / normalized_length * 100;
                                    p_val = comp_data.clusters{cluster}.P;
                                    
                                    fprintf('%-15s %-20s %-15s %-20.1f %-15.1f %-15.4f\n', ...
                                        '', comp_data.comparison, comp_data.test_type, ...
                                        start_pct, end_pct, p_val);
                                end
                            else
                                fprintf('%-15s %-20s %-15s %-20s %-15s %-15s\n', ...
                                    '', comp_data.comparison, comp_data.test_type, ...
                                    'None', 'signif.', 'p > 0.05');
                            end
                        end
                        fprintf('\n');
                    end
                end
            end
        end

end

end
    