function plotCombinedEMGPerSubjectWithSPM1DCycles(data_files)
    % Cette fonction trace les profils EMG moyens des 4 mouvements combinés 
    % pour tous les sujets, avec une figure distincte pour chaque muscle,
    % et réalise une analyse SPM1D pour comparer les groupes.
    % Cette fonction permet de sélectionner trois cycles représentatifs avec ginput 
    % qui sont ensuite moyennés en UN SEUL cycle moyen représentatif pour l'analyse SPM1D.
    %
    % Paramètres:
    %   data_files - Cellule contenant les chemins des fichiers .mat (jusqu'à 3)
    %               ou chaîne de caractères pour un seul fichier

    % Ajouter le chemin vers SPM1D
    addpath(genpath('C:\Users\Francalanci Hugo\Documents\MATLAB\Stage Sainte-Justine\HUG\Statistics\spm1dmatlab-master'))  

    % Vérifier le type de l'entrée et la convertir en cellule si nécessaire
    if ischar(data_files) || isstring(data_files)
        data_files = {data_files};
    end

    % Limiter à 3 fichiers maximum
    if length(data_files) > 3
        warning('Un maximum de 3 fichiers est supporté. Seuls les 3 premiers seront traités.');
        data_files = data_files(1:3);
    end

    % Nombre de fichiers à traiter
    num_files = length(data_files);

    % Définir les styles de ligne pour chaque fichier
    line_styles = {'-', '-', '-'};

    % Définir les couleurs pour chaque groupe (RGB)
    line_colors = cell(3, 1);
    fill_colors = cell(3, 1);

    % Couleurs par défaut pour les 3 groupes possibles
    % Noir pour Asympto, Rouge pour Preop, Bleu pour Postop
    default_colors = {
        [0, 0, 0],          % noir (Asympto)
        [1, 0, 0],          % rouge (Preop)
        [0, 0, 1]           % bleu (Postop)
    };

    % Couleurs pour les intervalles de confiance (plus claires)
    default_fill_colors = {
        [0.5, 0.5, 0.5],    % gris (Asympto)
        [1, 0.5, 0.5],      % rouge clair (Preop)
        [0.5, 0.5, 1]       % bleu clair (Postop)
    };

    % Pour stocker les données de chaque fichier
    all_data = cell(num_files, 1);
    file_labels = cell(num_files, 1);

    % Charger tous les fichiers
    for f = 1:num_files
        % Vérifier si le fichier existe
        if ~exist(data_files{f}, 'file')
            error('Le fichier de données spécifié n''existe pas: %s', data_files{f});
        end

        % Charger les données
        fprintf('Chargement des données depuis %s...\n', data_files{f});
        data = load(data_files{f});

        % Vérifier que la structure contient les données nécessaires
        if ~isfield(data, 'individual_data')
            error('Le fichier %s ne contient pas la structure "individual_data" attendue.', data_files{f});
        end

        % Stocker les données
        all_data{f} = data.individual_data;

        % Extraire le nom du fichier pour le label
        [~, file_name, ~] = fileparts(data_files{f});
        file_labels{f} = file_name;
    end

    % Déterminer les labels et couleurs pour les groupes
    % Par défaut: asympt, pre, post (si disponible)
    group_labels = cell(num_files, 1);

    % Initialiser les couleurs pour chaque fichier
    for f = 1:num_files
        line_colors{f} = default_colors{1};  % Défaut à noir
        fill_colors{f} = default_fill_colors{1};  % Défaut à gris
    end

    if num_files == 1
        group_labels{1} = 'Groupe';
    elseif num_files == 2
        if contains(lower(file_labels{1}), 'pre')
            group_labels{1} = 'Pre';
            line_colors{1} = default_colors{2};  % Rouge pour Pre
            fill_colors{1} = default_fill_colors{2};

            if contains(lower(file_labels{2}), 'post')
                group_labels{2} = 'Post';
                line_colors{2} = default_colors{3};  % Bleu pour Post
                fill_colors{2} = default_fill_colors{3};
            else
                group_labels{2} = 'Asympt';
                line_colors{2} = default_colors{1};  % Noir pour Asympt
                fill_colors{2} = default_fill_colors{1};
            end
        elseif contains(lower(file_labels{1}), 'asympt')
            group_labels{1} = 'Asympt';
            line_colors{1} = default_colors{1};  % Noir pour Asympt
            fill_colors{1} = default_fill_colors{1};

            if contains(lower(file_labels{2}), 'pre')
                group_labels{2} = 'Pre';
                line_colors{2} = default_colors{2};  % Rouge pour Pre
                fill_colors{2} = default_fill_colors{2};
            else
                group_labels{2} = 'Post';
                line_colors{2} = default_colors{3};  % Bleu pour Post
                fill_colors{2} = default_fill_colors{3};
            end
        else
            group_labels{1} = 'Groupe 1';
            group_labels{2} = 'Groupe 2';
        end
    elseif num_files == 3
        for f = 1:num_files
            if contains(lower(file_labels{f}), 'asympt')
                group_labels{f} = 'Asympt';
                line_colors{f} = default_colors{1};  % Noir pour Asympt
                fill_colors{f} = default_fill_colors{1};
            elseif contains(lower(file_labels{f}), 'pre')
                group_labels{f} = 'Pre';
                line_colors{f} = default_colors{2};  % Rouge pour Pre
                fill_colors{f} = default_fill_colors{2};
            elseif contains(lower(file_labels{f}), 'post')
                group_labels{f} = 'Post';
                line_colors{f} = default_colors{3};  % Bleu pour Post
                fill_colors{f} = default_fill_colors{3};
            else
                group_labels{f} = ['Groupe ' num2str(f)];
            end
        end
    end

    % Extraire les informations nécessaires du premier fichier pour initialisation
    time_normalized = all_data{1}.time;
    muscles_R = all_data{1}.muscles_R;
    muscles_L = all_data{1}.muscles_L;

    % Nombre de muscles
    nb_muscles = length(muscles_R);
    if nb_muscles ~= length(muscles_L)
        warning('Le nombre de muscles diffère entre les côtés droit et gauche.');
        nb_muscles = min(length(muscles_R), length(muscles_L));
    end

    % Créer une structure pour stocker les données de muscle par fichier pour l'analyse SPM1D
    all_muscle_data = cell(nb_muscles, 1);
    for m = 1:nb_muscles
        all_muscle_data{m} = cell(num_files, 1);
    end

    % Pour chaque muscle, créer une figure séparée pour les profils EMG
    for m = 1:nb_muscles
        % Créer une nouvelle figure pour ce muscle
        figure('Name', sprintf('Muscle %d: %s', m, muscles_R{m}), 'Color', 'white', 'Position', [100+m*50, 100+m*30, 800, 500]);
        hold on;

        % Légende pour cette figure
        legend_handles = [];
        legend_labels = {};

        % Traiter chaque fichier
        for f = 1:num_files
            % Extraire les données de ce fichier
            data = all_data{f};

            % Extraire les informations spécifiques
            functional_labels = data.functional_labels;
            subject_data_R = data.subject_data_R;
            subject_data_L = data.subject_data_L;
            subject_ids_R = data.subject_ids_R;
            subject_ids_L = data.subject_ids_L;

            % Nombre de mouvements fonctionnels
            nb_functionals = length(functional_labels);

            % Nombre de sujets pour chaque côté
            nb_subjects_R = length(subject_ids_R);
            nb_subjects_L = length(subject_ids_L);

            % Matrices pour stocker toutes les moyennes de sujets
            all_means_R = zeros(nb_subjects_R, length(time_normalized));
            all_means_L = zeros(nb_subjects_L, length(time_normalized));

            % Calculer la moyenne des 4 mouvements pour chaque sujet côté droit
            for s = 1:nb_subjects_R
                subject_mean = zeros(1, length(time_normalized));

                for func = 1:nb_functionals
                    % Extraire les données de ce sujet pour ce mouvement
                    subj_data = subject_data_R{s, func};

                    % Vérifier les dimensions et extraire les données du muscle
                    if size(subj_data, 1) == nb_muscles
                        % Les muscles sont en ligne
                        muscle_data = subj_data(m, :);
                    else
                        % Les muscles sont en colonne
                        muscle_data = subj_data(:, m)';
                    end

                    % Ajouter à la moyenne
                    subject_mean = subject_mean + muscle_data;
                end

                % Calculer la moyenne des mouvements pour ce sujet
                subject_mean = subject_mean / nb_functionals;

                % Stocker dans le tableau
                all_means_R(s, :) = subject_mean;
            end

            % Calculer la moyenne des 4 mouvements pour chaque sujet côté gauche
            for s = 1:nb_subjects_L
                subject_mean = zeros(1, length(time_normalized));

                for func = 1:nb_functionals
                    % Extraire les données de ce sujet pour ce mouvement
                    subj_data = subject_data_L{s, func};

                    % Vérifier les dimensions et extraire les données du muscle
                    if size(subj_data, 1) == nb_muscles
                        % Les muscles sont en ligne
                        muscle_data = subj_data(m, :);
                    else
                        % Les muscles sont en colonne
                        muscle_data = subj_data(:, m)';
                    end

                    % Ajouter à la moyenne
                    subject_mean = subject_mean + muscle_data;
                end

                % Calculer la moyenne des mouvements pour ce sujet
                subject_mean = subject_mean / nb_functionals;

                % Stocker dans le tableau
                all_means_L(s, :) = subject_mean;
            end

            % Combiner les données des côtés droit et gauche
            all_means_combined = [all_means_R; all_means_L];

            % Stocker les données pour l'analyse SPM1D
            all_muscle_data{m}{f} = all_means_combined;

            % Calculer la moyenne finale et l'écart-type de tous les sujets
            mean_data = mean(all_means_combined, 1, 'omitnan');
            std_data = std(all_means_combined, 0, 1, 'omitnan');

            % Calculer l'intervalle de confiance à 95%
            n_subjects = size(all_means_combined, 1);
            sem = std_data / sqrt(n_subjects);
            t_crit = tinv(0.975, n_subjects - 1);  % 95% de confiance
            ci_lower = mean_data - t_crit * sem;
            ci_upper = mean_data + t_crit * sem;

            % Tracer la ligne moyenne avec le style approprié
            h_line = plot(time_normalized, mean_data, line_styles{f}, 'LineWidth', 2.5, 'Color', line_colors{f});

            % Tracer l'intervalle de confiance avec fill
            x_fill = [time_normalized, fliplr(time_normalized)];
            y_fill = [ci_upper, fliplr(ci_lower)];
            fill(x_fill, y_fill, fill_colors{f}, 'EdgeColor', 'none', 'FaceAlpha', 0.2);

            % Ajouter à la légende
            legend_handles = [legend_handles, h_line];

            % Créer une étiquette avec le nom du fichier et le nombre de sujets
            legend_label = sprintf('%s (n=%d)', group_labels{f}, n_subjects);
            legend_labels{end+1} = legend_label;
        end

        % Définir le titre de la figure
        if isequal(muscles_R{m}, muscles_L{m})
            title_str = muscles_R{m};
        else
            title_str = [muscles_R{m} ' / ' muscles_L{m}];
        end
        title(title_str, 'FontWeight', 'bold', 'FontSize', 14);

        % Axes et grille
        xlabel('Temps normalisé (%)', 'FontSize', 12);
        ylabel('EMG normalisé (% MVC)', 'FontSize', 12);
        xlim([min(time_normalized), max(time_normalized)]);
        grid off;
        box on;

        % Ajouter la légende
        if ~isempty(legend_handles)
            legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 12);
        end

        % Ajuster la mise en page
        set(gca, 'FontSize', 12);
    end

   % ========== SYSTÈME DE SAUVEGARDE/CHARGEMENT DES SÉLECTIONS ==========

% Créer un nom de fichier de sauvegarde basé sur les fichiers d'entrée
save_filename = 'cycle_selections';
for f = 1:num_files
    [~, name, ~] = fileparts(data_files{f});
    save_filename = [save_filename '_' name];
end
save_filename = [save_filename '.mat'];

% Variables pour stocker les sélections
saved_selections = struct();
load_previous = false;

% Vérifier s'il existe déjà un fichier de sélections
if exist(save_filename, 'file')
    fprintf('\n=== FICHIER DE SÉLECTIONS DÉTECTÉ ===\n');
    fprintf('Un fichier de sélections existe déjà : %s\n', save_filename);
    
    % Charger les sélections précédentes
    try
        loaded_data = load(save_filename);
        if isfield(loaded_data, 'cycle_selections')
            saved_selections = loaded_data.cycle_selections;
            
            % Vérifier la compatibilité
            if isfield(saved_selections, 'nb_muscles') && ...
               isfield(saved_selections, 'num_files') && ...
               isfield(saved_selections, 'muscle_names') && ...
               saved_selections.nb_muscles == nb_muscles && ...
               saved_selections.num_files == num_files
                
                % Vérifier que les noms de muscles correspondent
                muscles_match = true;
                for m = 1:nb_muscles
                    if ~strcmp(saved_selections.muscle_names{m}, muscles_R{m})
                        muscles_match = false;
                        break;
                    end
                end
                
                if muscles_match
                    fprintf('✓ Sélections compatibles trouvées !\n');
                    fprintf('  - %d muscles\n', saved_selections.nb_muscles);
                    fprintf('  - %d fichiers\n', saved_selections.num_files);
                    fprintf('  - Date de création : %s\n', saved_selections.creation_date);
                    
                    % Demander à l'utilisateur
                    user_choice = input('Voulez-vous utiliser les sélections précédentes ? (y/n): ', 's');
                    if strcmpi(user_choice, 'y') || strcmpi(user_choice, 'yes')
                        load_previous = true;
                        fprintf('→ Chargement des sélections précédentes...\n');
                    else
                        fprintf('→ Nouvelle sélection manuelle...\n');
                    end
                else
                    fprintf('⚠ Les noms de muscles ne correspondent pas. Nouvelle sélection requise.\n');
                end
            else
                fprintf('⚠ Structure incompatible. Nouvelle sélection requise.\n');
            end
        end
    catch ME
        fprintf('⚠ Erreur lors du chargement : %s\n', ME.message);
        fprintf('→ Nouvelle sélection manuelle...\n');
    end
    fprintf('=====================================\n\n');
end

% Structure pour stocker les cycles sélectionnés et les indices
selected_cycles = cell(nb_muscles, num_files);
cycle_indices = cell(nb_muscles, num_files);  % Nouveauté : stocker les indices
representative_cycles = cell(nb_muscles, num_files);
normalized_length = 100;

% Boucle de sélection des cycles 
for m = 1:nb_muscles
    fprintf('\n===== Sélection des cycles pour le muscle %d: %s =====\n', m, muscles_R{m});

    for f = 1:num_files
        % Vérifier s'il faut charger ou sélectionner
        if load_previous && isfield(saved_selections, 'selections') && ...
           size(saved_selections.selections, 1) >= m && ...
           size(saved_selections.selections, 2) >= f && ...
           ~isempty(saved_selections.selections{m, f})
            
            % CHARGER LES SÉLECTIONS PRÉCÉDENTES
            fprintf('Chargement de la sélection précédente pour %s - %s...\n', muscles_R{m}, group_labels{f});
            
            % Récupérer les indices sauvegardés
            saved_indices = saved_selections.selections{m, f};
            cycle_indices{m, f} = saved_indices;
            
            % Recréer les cycles à partir des indices
            cycles = cell(3, 1);
            for c = 1:3
                start_idx = saved_indices(c, 1);
                end_idx = saved_indices(c, 2);
                
                if start_idx > 0 && end_idx > 0 && start_idx <= end_idx && ...
                   end_idx <= length(time_normalized)
                    
                    cycle_range = start_idx:end_idx;
                    cycles{c} = all_muscle_data{m}{f}(:, cycle_range);
                else
                    fprintf('⚠ Indices invalides pour le cycle %d, utilisation du signal complet\n', c);
                    cycles{c} = all_muscle_data{m}{f};
                end
            end
            
            selected_cycles{m, f} = cycles;
            
            % Afficher un graphique de confirmation (optionnel)
            figure('Name', sprintf('Confirmation - Muscle %d: %s (%s)', m, muscles_R{m}, group_labels{f}), ...
                  'Color', 'white', 'Position', [200+f*50, 200+f*30, 800, 500]);
            
            mean_data = mean(all_muscle_data{m}{f}, 1, 'omitnan');
            plot(time_normalized, mean_data, line_styles{f}, 'LineWidth', 2, 'Color', line_colors{f});
            hold on;
            
            % Afficher les cycles sélectionnés
            for c = 1:3
                if ~isempty(cycles{c}) && size(cycles{c}, 2) > 1
                    start_idx = saved_indices(c, 1);
                    end_idx = saved_indices(c, 2);
                    cycle_range = start_idx:end_idx;
                    plot(time_normalized(cycle_range), mean(cycles{c}, 1, 'omitnan'), 'LineWidth', 2, 'Color', [0.2, 0.6, 0.2]);
                end
            end
            
            title(sprintf('Muscle %s - %s: Cycles chargés', muscles_R{m}, group_labels{f}), 'FontSize', 14);
            xlabel('Temps normalisé (%)', 'FontSize', 12);
            ylabel('EMG normalisé (% MVC)', 'FontSize', 12);
            legend('Signal complet', 'Cycles chargés', 'Location', 'northeast');
            
        else
            % SÉLECTION MANUELLE 
            figure('Name', sprintf('Sélection des cycles - Muscle %d: %s (%s)', m, muscles_R{m}, group_labels{f}), ...
                  'Color', 'white', 'Position', [200+f*50, 200+f*30, 800, 500]);

            mean_data = mean(all_muscle_data{m}{f}, 1, 'omitnan');
            plot(time_normalized, mean_data, line_styles{f}, 'LineWidth', 2, 'Color', line_colors{f});

            title(sprintf('Muscle %s - %s: Sélectionnez 3 cycles représentatifs', muscles_R{m}, group_labels{f}), 'FontSize', 14);
            xlabel('Temps normalisé (%)', 'FontSize', 12);
            ylabel('EMG normalisé (% MVC)', 'FontSize', 12);
            xlim([min(time_normalized), max(time_normalized)]);
            grid off;

            annotation('textbox', [0.15, 0.01, 0.7, 0.05], 'String', ...
                'Cliquez pour sélectionner le début et la fin de chaque cycle (6 points au total)', ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

            fprintf('Veuillez sélectionner les points de début et de fin de 3 cycles pour %s - %s\n', muscles_R{m}, group_labels{f});
            [x_points, ~] = ginput(6);
            x_points = sort(x_points);

            if length(x_points) ~= 6
                warning('Nombre incorrect de points sélectionnés. Utilisation de l''ensemble du signal.');
                selected_cycles{m, f} = all_muscle_data{m}{f};
                cycle_indices{m, f} = [];  % Pas d'indices valides
                continue;
            end

            % Extraire les cycles et stocker les indices
            cycles = cell(3, 1);
            indices_for_save = zeros(3, 2);  % [start, end] pour chaque cycle
            hold on;

            for c = 1:3
                start_idx = find(time_normalized >= x_points(c*2-1), 1, 'first');
                end_idx = find(time_normalized <= x_points(c*2), 1, 'last');

                if isempty(start_idx) || isempty(end_idx) || start_idx >= end_idx
                    warning('Points de cycle invalides. Utilisation de l''ensemble du signal.');
                    selected_cycles{m, f} = all_muscle_data{m}{f};
                    cycle_indices{m, f} = [];
                    continue;
                end

                % Stocker les indices
                indices_for_save(c, 1) = start_idx;
                indices_for_save(c, 2) = end_idx;

                cycle_range = start_idx:end_idx;
                cycles{c} = all_muscle_data{m}{f}(:, cycle_range);

                plot(time_normalized(cycle_range), mean(cycles{c}, 1, 'omitnan'), 'LineWidth', 2, 'Color', [0.2, 0.6, 0.2]);
            end

            selected_cycles{m, f} = cycles;
            cycle_indices{m, f} = indices_for_save;

            title(sprintf('Muscle %s - %s: 3 cycles sélectionnés', muscles_R{m}, group_labels{f}), 'FontSize', 14);
            legend('Signal complet', 'Cycles sélectionnés', 'Location', 'northeast');
        end
        
        % Créer le cycle moyen représentatif 
        num_subjects = size(all_muscle_data{m}{f}, 1);
        single_representative_cycle = zeros(num_subjects, normalized_length);
        
        cycles = selected_cycles{m, f};  % Utiliser les cycles (chargés ou sélectionnés)
        
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

% SAUVEGARDER LES NOUVELLES SÉLECTIONS
fprintf('\n=== SAUVEGARDE DES SÉLECTIONS ===\n');
try
    % Préparer la structure de sauvegarde
    selections_to_save = struct();
    selections_to_save.nb_muscles = nb_muscles;
    selections_to_save.num_files = num_files;
    selections_to_save.muscle_names = muscles_R;
    selections_to_save.group_labels = group_labels;
    selections_to_save.file_labels = file_labels;
    selections_to_save.creation_date = datestr(now);
    selections_to_save.selections = cycle_indices;  % Stocker les indices
    
    % Sauvegarder
    cycle_selections = selections_to_save;  % Variable pour le fichier .mat
    save(save_filename, 'cycle_selections');
    
    fprintf('✓ Sélections sauvegardées dans : %s\n', save_filename);
    fprintf('  - %d muscles x %d fichiers\n', nb_muscles, num_files);
    fprintf('  - Date : %s\n', datestr(now));
catch ME
    fprintf('⚠ Erreur lors de la sauvegarde : %s\n', ME.message);
end
fprintf('=================================\n\n');

    % ========== NOUVELLE FIGURE POUR ARTICLE ==========
    if num_files > 1
        % Créer la figure finale pour l'article
        figure('Name', 'Figure Article - Cycles représentatifs et analyses SPM1D', ...
               'Color', 'white', 'Position', [50, 50, 1500, 800]);
        
        % Structure pour stocker les résultats statistiques
        spm_results = struct();
        
        % Temps normalisé pour le cycle représentatif (0 à 100%)
        cycle_time_percent = linspace(0, 100, normalized_length);
        
        % Pour chaque muscle, créer les subplots
        for m = 1:nb_muscles
            % Déterminer le nom du muscle
            if isequal(muscles_R{m}, muscles_L{m})
                muscle_name = muscles_R{m};
            else
                muscle_name = [muscles_R{m} '/' muscles_L{m}];
            end
            
            % Subplot pour les cycles EMG 
            subplot(2, 3, m);
            hold on;
            
            % Tracer les cycles pour chaque condition
            legend_handles = [];
            legend_labels = {};
            
            for f = 1:num_files
                if ~isempty(representative_cycles{m, f})
                    % Calculer la moyenne et l'intervalle de confiance
                    mean_data = mean(representative_cycles{m, f}, 1, 'omitnan');
                    std_data = std(representative_cycles{m, f}, 0, 1, 'omitnan');
                    n_subjects = size(representative_cycles{m, f}, 1);
                    
                    sem = std_data / sqrt(n_subjects);
                    t_crit = tinv(0.975, n_subjects - 1);
                    ci_lower = mean_data - t_crit * sem;
                    ci_upper = mean_data + t_crit * sem;
                    
                    % Tracer la ligne moyenne
                    h_line = plot(cycle_time_percent, mean_data, line_styles{f}, ...
                        'LineWidth', 2.5, 'Color', line_colors{f});
                    
                    % Tracer l'intervalle de confiance
                    x_fill = [cycle_time_percent, fliplr(cycle_time_percent)];
                    y_fill = [ci_upper, fliplr(ci_lower)];
                    fill(x_fill, y_fill, fill_colors{f}, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
                    
                    % Ajouter à la légende
                    legend_handles = [legend_handles, h_line];
                    legend_labels{end+1} = sprintf('%s (n=%d)', group_labels{f}, n_subjects);
                end
            end
            
            % Configuration du subplot EMG
            title(muscle_name, 'FontWeight', 'bold', 'FontSize', 16);
            xlabel('Cycle (%)', 'FontSize', 14);
            if m == 1 || m == 4
               ylabel('Normalized EMG (%)', 'FontSize', 14);
            else
               set(gca, 'YTickLabel', []); % Supprime les labels des ticks Y
               ylabel(''); % Supprime le label de l'axe Y
            end
            xlim([0, 100]);
            ylim([0, 60]);  % Limites Y fixes pour tous les graphiques
            grid off;
            box on;

            if num_files >= 2
            % Définir les valeurs pour chaque groupe
            asympt_line = 58.5;
            pre_line = 54.5;
            post_line = 49.5;
    
    % Trouver quels groupes sont présents
    has_asympt = any(contains(group_labels, 'Asympt'));
    has_pre = any(contains(group_labels, 'Pre'));
    has_post = any(contains(group_labels, 'Post'));
    
    % Tracer les lignes selon les groupes présents
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

            
            % Légende seulement sur le premier subplot
            if m == 1 && ~isempty(legend_handles)
                legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 13);
            end

% ========== ANALYSE SPM1D - BARRES SOUS LE GRAPHIQUE ==========

% Initialiser les résultats pour ce muscle
spm_muscle_results = struct();

% Effectuer les comparaisons par paires
comparisons = {};
bar_colors = {};

% ========== ANOVA À UN FACTEUR - SECTION CORRIGÉE ==========
% Préparer les données pour l'ANOVA
anova_data = {};
anova_labels = {};

for f = 1:num_files
    if ~isempty(representative_cycles{m, f})
        anova_data{end+1} = representative_cycles{m, f};
        anova_labels{end+1} = group_labels{f};
    end
end

% Effectuer l'ANOVA SPM1D si au moins 2 groupes
anova_result = struct();
if length(anova_data) >= 2
    try
        % Vérifier les dimensions des données
        fprintf('Dimensions des données pour ANOVA - Muscle %s:\n', muscle_name);
        for i = 1:length(anova_data)
            fprintf('  Groupe %d (%s): %d sujets x %d points temporels\n', ...
                i, anova_labels{i}, size(anova_data{i}, 1), size(anova_data{i}, 2));
        end
        
        % Méthode 1: Essayer avec les données sous forme de cellule (format standard SPM1D)
        try
            spm_anova = spm1d.stats.anova1(anova_data);
            fprintf('ANOVA réussie avec format cellule\n');
        catch
            % Méthode 2: Convertir en format matriciel avec vecteur de groupes
            fprintf('Format cellule échoué, essai avec format matriciel...\n');
            
            % Combiner toutes les données en une seule matrice
            all_subjects_data = [];
            group_vector = [];
            
            for i = 1:length(anova_data)
                all_subjects_data = [all_subjects_data; anova_data{i}];
                group_vector = [group_vector; repmat(i, size(anova_data{i}, 1), 1)];
            end
            
            % Appel ANOVA avec matrice et vecteur de groupes
            spm_anova = spm1d.stats.anova1(all_subjects_data, group_vector);
            fprintf('ANOVA réussie avec format matriciel\n');
        end
        
        % Inférence statistique
        spmi_anova = spm_anova.inference(0.05, 'interp', true);
        
        % Stocker les résultats de l'ANOVA
        anova_result.F_stat = spm_anova.z; % Statistique F
        anova_result.p_values = spmi_anova.p;
        anova_result.significant_clusters = spmi_anova.clusters;
        anova_result.global_significant = ~isempty(spmi_anova.clusters);
        
        fprintf('ANOVA - Muscle %s: ', muscle_name);
        if anova_result.global_significant
            fprintf('Effet significatif détecté (p < 0.05)\n');
            fprintf('  Nombre de clusters significatifs: %d\n', length(spmi_anova.clusters));
        else
            fprintf('Aucun effet significatif (p > 0.05)\n');
        end
        
    catch ME
        fprintf('Erreur lors de l''ANOVA SPM1D pour %s: %s\n', muscle_name, ME.message);
        fprintf('Détails de l''erreur: %s\n', ME.getReport);
        anova_result.global_significant = false;
        
        % Diagnostic supplémentaire
        fprintf('Diagnostic des données:\n');
        for i = 1:length(anova_data)
            if ~isempty(anova_data{i})
                fprintf('  Groupe %d: Taille = [%d x %d], Min = %.3f, Max = %.3f, NaN = %d\n', ...
                    i, size(anova_data{i}, 1), size(anova_data{i}, 2), ...
                    min(anova_data{i}(:), [], 'omitnan'), max(anova_data{i}(:), [], 'omitnan'), ...
                    sum(isnan(anova_data{i}(:))));
            else
                fprintf('  Groupe %d: Données vides\n', i);
            end
        end
    end
else
    fprintf('ANOVA impossible - Moins de 2 groupes disponibles pour %s\n', muscle_name);
    anova_result.global_significant = false;
end

if num_files == 3
   % Toutes les comparaisons possibles
   comparisons = {'Asympt vs Pre', 'Asympt vs Post', 'Pre vs Post'};

   % Trouver les indices correspondants
   asympt_idx = find(contains(group_labels, 'Asympt'));
   pre_idx = find(contains(group_labels, 'Pre'));
   post_idx = find(contains(group_labels, 'Post'));

   pairs = [asympt_idx pre_idx; asympt_idx post_idx; pre_idx post_idx];
   bar_colors = {[1, 0, 0], [0, 0, 1], [1, 0.5, 0]}; % rouge, bleu, orange
    
elseif num_files == 2
    % Une seule comparaison
    comparisons = {sprintf('%s vs %s', group_labels{1}, group_labels{2})};
    pairs = [1 2];
    
    % Déterminer la couleur selon le type de comparaison
    if contains(group_labels{1}, 'Asympt') || contains(group_labels{2}, 'Asympt')
        bar_colors = {[0, 0, 0]}; % noir
    elseif contains(group_labels{1}, 'Pre') || contains(group_labels{2}, 'Pre')
        bar_colors = {[1, 0, 0]}; % rouge
    else
        bar_colors = {[0, 0, 1]}; % bleu
    end
end

% Variables pour stocker les p-values
all_p_values = {};

% Obtenir les limites actuelles du graphique
current_ylim = ylim;
y_min = 0;
y_max = 60;

% ========== TESTS T CONDITIONNELS ==========
% Effectuer les tests t seulement si l'ANOVA est significative
if anova_result.global_significant
    fprintf('ANOVA significative - Procédure aux tests post-hoc pour %s...\n', muscle_name);
    
    % Effectuer chaque comparaison
    for comp = 1:size(pairs, 1)
        i = pairs(comp, 1);
        j = pairs(comp, 2);
    
    if ~isempty(representative_cycles{m, i}) && ~isempty(representative_cycles{m, j})
        % Données alignées
        data1 = representative_cycles{m, i};
        data2 = representative_cycles{m, j};
        
        % Test statistique
        try
            if size(data1, 1) == size(data2, 1) && ...
               (contains(comparisons{comp}, 'Pre vs Post'))
                % Test apparié pour Pre vs Post
                spm = spm1d.stats.ttest_paired(data1, data2);
                test_type = 'paired';
            else
                % Test indépendant
                spm = spm1d.stats.ttest2(data1, data2);
                test_type = 'independent';
            end
            
            % Inférence statistique
            spmi = spm.inference(0.05, 'two_tailed', true, 'interp', true);
            
            % Extraire les régions significatives
            if ~isempty(spmi.clusters)
                for cluster = 1:length(spmi.clusters)
                    start_pct = (spmi.clusters{cluster}.endpoints(1) - 1) / normalized_length * 100;
                    end_pct = (spmi.clusters{cluster}.endpoints(2) - 1) / normalized_length * 100;
                    
                    % Position des barres sous le graphique principal
                    bar_height = 1; % Hauteur fixe
                    y_position = 5 - (comp * 1); % positions
                    
                    % Dessiner la barre horizontale sous le graphique
                    rectangle('Position', [start_pct, y_position, end_pct - start_pct, bar_height], ...
                              'FaceColor', bar_colors{comp}, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                end
                
                % Stocker les p-values (même code qu'avant)
                cluster_info = struct();
                cluster_info.comparison = comparisons{comp};
                cluster_info.test_type = test_type;
                cluster_info.clusters = spmi.clusters;
                cluster_info.p_values = spmi.p;
                
                all_p_values{end+1} = cluster_info;
            end
            
        catch ME
            warning('Erreur lors de l''analyse SPM1D pour %s - %s: %s', ...
                muscle_name, comparisons{comp}, ME.message);
        end
    end
    end

    else
    fprintf('ANOVA non significative - Tests post-hoc non nécessaires pour %s\n', muscle_name);
    all_p_values = {}; % Aucun test t effectué
end

% Stocker les résultats pour ce muscle
spm_muscle_results.muscle_name = muscle_name;
spm_muscle_results.anova_result = anova_result;
spm_muscle_results.comparisons = all_p_values;
spm_results.(sprintf('muscle_%d', m)) = spm_muscle_results;

% Ajouter une légende pour les barres sur le dernier subplot
if m == nb_muscles && ~isempty(comparisons)
    % Créer des éléments invisibles pour la légende des barres
    legend_handles_spm = [];
    for comp = 1:length(comparisons)
        h_dummy = plot(NaN, NaN, 's', 'MarkerFaceColor', bar_colors{comp}, ...
               'MarkerEdgeColor', 'none', 'MarkerSize', 8);
        legend_handles_spm = [legend_handles_spm, h_dummy];
    end
    
    % Ajouter la légende des barres
    legend_combined = [legend_handles, legend_handles_spm];
    legend_labels_combined = [legend_labels, comparisons];
    legend(legend_combined, legend_labels_combined, 'Location', 'best', 'FontSize', 13);
end
        
        % Titre global
        sgtitle('Cycles représentatifs EMG et Analyses statistiques SPM1D', 'FontSize', 16, 'FontWeight', 'bold');
            
        % ========== TABLEAU DES P-VALUES ==========
        fprintf('\n========== RÉSULTATS DES ANALYSES SPM1D ==========\n');
        fprintf('Tableau des p-values et régions significatives:\n\n');
        
                % En-tête du tableau
                fprintf('%-15s %-20s %-15s %-20s %-15s %-15s\n', ...
                    'Muscle', 'Comparaison', 'Type Test', 'Début Cycle (%)', 'Fin Cycle (%)', 'p-value');
                fprintf(repmat('-', 1, 120));
                fprintf('\n');
                
                % Afficher les résultats pour chaque muscle
                for m = 1:nb_muscles
                    muscle_field = sprintf('muscle_%d', m);
                    if isfield(spm_results, muscle_field)
                        muscle_results = spm_results.(muscle_field);
                        
                        if ~isempty(muscle_results.comparisons)
                            % Afficher le nom du muscle
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
                                        'Aucune', 'signif.', 'p > 0.05');
                                end
                            end
                            fprintf('\n');
                        end
                    end
                end
                
                fprintf('========== FIN DES RÉSULTATS SPM1D ==========\n\n');
                
% ========== SAUVEGARDE EXCEL ==========
% Créer un tableau récapitulatif pour Excel
excel_data = {};
row_idx = 1;

% En-têtes
excel_data{row_idx, 1} = 'Muscle';
excel_data{row_idx, 2} = 'ANOVA_Significative';
excel_data{row_idx, 3} = 'ANOVA_P_min';
excel_data{row_idx, 4} = 'ANOVA_P_global'; 
excel_data{row_idx, 5} = 'Type_Test';      
excel_data{row_idx, 6} = 'Comparaison';
excel_data{row_idx, 7} = 'Debut_Cycle_Pourcent';
excel_data{row_idx, 8} = 'Fin_Cycle_Pourcent';
excel_data{row_idx, 9} = 'P_Value';
excel_data{row_idx, 10} = 'Significatif';
row_idx = row_idx + 1;

% Remplir les données pour chaque muscle
for m = 1:nb_muscles
    muscle_field = sprintf('muscle_%d', m);
    if isfield(spm_results, muscle_field)
        muscle_results = spm_results.(muscle_field);

        % Informations ANOVA
        muscle_anova = muscle_results.anova_result;
        anova_significant = muscle_anova.global_significant;

        % P-value minimale de l'ANOVA
if anova_significant && ~isempty(muscle_anova.significant_clusters)
    p_values = [];
    for k = 1:length(muscle_anova.significant_clusters)
        p_values = [p_values, muscle_anova.significant_clusters{k}.P];
    end
    anova_p_min = min(p_values);
    anova_p_global = min(muscle_anova.p_values);  % <-- NOUVEAU
else
    anova_p_min = NaN;
    anova_p_global = NaN;  % <-- NOUVEAU
end

        % Si pas de tests t (ANOVA non significative)
        if isempty(muscle_results.comparisons)
            excel_data{row_idx, 1} = muscle_results.muscle_name;
            excel_data{row_idx, 2} = anova_significant;
            excel_data{row_idx, 3} = anova_p_min;
            excel_data{row_idx, 4} = anova_p_global;  
            excel_data{row_idx, 5} = 'ANOVA_seulement';  
            excel_data{row_idx, 6} = 'Aucun_test_posthoc';
            excel_data{row_idx, 7} = NaN;
            excel_data{row_idx, 8} = NaN;
            excel_data{row_idx, 9} = NaN;
            excel_data{row_idx, 10} = false;
            row_idx = row_idx + 1;
        else
            % Tests t effectués
            for comp_idx = 1:length(muscle_results.comparisons)
                comp_data = muscle_results.comparisons{comp_idx};

                if ~isempty(comp_data.clusters)
                    % Clusters significatifs
                    for cluster = 1:length(comp_data.clusters)
                        start_pct = (comp_data.clusters{cluster}.endpoints(1) - 1) / normalized_length * 100;
                        end_pct = (comp_data.clusters{cluster}.endpoints(2) - 1) / normalized_length * 100;
                        p_val = comp_data.clusters{cluster}.P;

excel_data{row_idx, 1} = muscle_results.muscle_name;
excel_data{row_idx, 2} = anova_significant;
excel_data{row_idx, 3} = anova_p_min;
excel_data{row_idx, 4} = anova_p_global;  
excel_data{row_idx, 5} = comp_data.test_type;
excel_data{row_idx, 6} = comp_data.comparison;
excel_data{row_idx, 7} = start_pct;
excel_data{row_idx, 8} = end_pct;
excel_data{row_idx, 9} = p_val;
excel_data{row_idx, 10} = true;
                        row_idx = row_idx + 1;
                    end
                else
                    % Pas de clusters significatifs
excel_data{row_idx, 1} = muscle_results.muscle_name;
excel_data{row_idx, 2} = anova_significant;
excel_data{row_idx, 3} = anova_p_min;
excel_data{row_idx, 4} = anova_p_global;  
excel_data{row_idx, 5} = comp_data.test_type;
excel_data{row_idx, 6} = comp_data.comparison;
excel_data{row_idx, 7} = NaN;
excel_data{row_idx, 8} = NaN;
excel_data{row_idx, 9} = NaN;
excel_data{row_idx, 10} = false;
                    row_idx = row_idx + 1;
                end
            end
        end
    end
end

% Sauvegarder en Excel
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
excel_filename = sprintf('SPM1D_Analysis_Results_%s.xlsx', timestamp);

% Convertir en table pour Excel
T = cell2table(excel_data(2:end, :), 'VariableNames', excel_data(1, :));
writetable(T, excel_filename, 'Sheet', 'Resultats_SPM1D');

fprintf('Résultats sauvegardés dans: %s\n', excel_filename);

    end

    % ========== EXTRACTION DES VALEURS EMG POUR JASP ==========
% Créer un tableau avec toutes les valeurs EMG des cycles moyens
fprintf('\nExtraction des valeurs EMG des cycles moyens pour analyse JASP...\n');

% Préparer la structure de données pour JASP
jasp_data = {};
col_idx = 1;

% En-têtes : Sujet, Groupe, puis chaque muscle à chaque point temporel
jasp_data{1, col_idx} = 'Sujet_ID'; col_idx = col_idx + 1;
jasp_data{1, col_idx} = 'Groupe'; col_idx = col_idx + 1;

% Créer les en-têtes pour chaque muscle et chaque point temporel
for m = 1:nb_muscles
    muscle_name = muscles_R{m};
    % Simplifier le nom du muscle pour Excel (enlever les caractères spéciaux)
    muscle_name = regexprep(muscle_name, '[^A-Za-z0-9]', '_');

    for t = 1:normalized_length
        jasp_data{1, col_idx} = sprintf('%s_T%d', muscle_name, t);
        col_idx = col_idx + 1;
    end
end

% Remplir les données ligne par ligne
row_idx = 2;
sujet_global_id = 1;

for f = 1:num_files
    group_name = group_labels{f};

    % Déterminer le nombre de sujets pour ce groupe
    if ~isempty(representative_cycles{1, f})
        n_subjects_group = size(representative_cycles{1, f}, 1);

        for s = 1:n_subjects_group
            col_idx = 1;

            % ID du sujet et groupe
            jasp_data{row_idx, col_idx} = sujet_global_id; col_idx = col_idx + 1;
            jasp_data{row_idx, col_idx} = group_name; col_idx = col_idx + 1;

            % Données EMG pour chaque muscle
            for m = 1:nb_muscles
                if ~isempty(representative_cycles{m, f})
                    % Vérifier que le sujet existe dans les données de ce muscle
                    if size(representative_cycles{m, f}, 1) >= s
                        muscle_data = representative_cycles{m, f}(s, :);

                        % Ajouter chaque point temporel
                        for t = 1:normalized_length
                            jasp_data{row_idx, col_idx} = muscle_data(t);
                            col_idx = col_idx + 1;
                        end
                    else
                        % Si pas de données pour ce sujet/muscle, remplir avec NaN
                        for t = 1:normalized_length
                            jasp_data{row_idx, col_idx} = NaN;
                            col_idx = col_idx + 1;
                        end
                    end
                else
                    % Si pas de données pour ce muscle, remplir avec NaN
                    for t = 1:normalized_length
                        jasp_data{row_idx, col_idx} = NaN;
                        col_idx = col_idx + 1;
                    end
                end
            end

            row_idx = row_idx + 1;
            sujet_global_id = sujet_global_id + 1;
        end
    end
end

% Sauvegarder le tableau JASP
jasp_filename = sprintf('EMG_Cycles_Moyens_JASP_%s.xlsx', timestamp);

% Convertir en table pour Excel
T_jasp = cell2table(jasp_data(2:end, :), 'VariableNames', jasp_data(1, :));
writetable(T_jasp, jasp_filename, 'Sheet', 'Donnees_EMG_Cycles');

fprintf('Données EMG des cycles moyens sauvegardées dans: %s\n', jasp_filename);
fprintf('Format : %d sujets x %d variables (%d muscles x %d points temporels + identificateurs)\n', ...
    size(T_jasp, 1), size(T_jasp, 2), nb_muscles, normalized_length);

% Préparation
jasp_data_mean = {};
jasp_data_mean{1,1} = 'Sujet_ID';
jasp_data_mean{1,2} = 'Groupe';
for m = 1:nb_muscles
    muscle_name = regexprep(muscles_R{m}, '[^A-Za-z0-9]', '_');
    jasp_data_mean{1, m+2} = muscle_name;
end

row_idx = 2;
sujet_global_id = 1;

for f = 1:num_files
    group_name = group_labels{f};
    if ~isempty(representative_cycles{1,f})
        n_subjects = size(representative_cycles{1,f}, 1);
        for s = 1:n_subjects
            jasp_data_mean{row_idx,1} = sujet_global_id;
            jasp_data_mean{row_idx,2} = group_name;
            for m = 1:nb_muscles
                if ~isempty(representative_cycles{m,f}) && size(representative_cycles{m,f},1) >= s
                    emg_values = representative_cycles{m,f}(s,:);
                    jasp_data_mean{row_idx, m+2} = mean(emg_values, 'omitnan');
                else
                    jasp_data_mean{row_idx, m+2} = NaN;
                end
            end
            row_idx = row_idx + 1;
            sujet_global_id = sujet_global_id + 1;
        end
    end
end

% Sauvegarde
T_jasp_mean = cell2table(jasp_data_mean(2:end,:), 'VariableNames', jasp_data_mean(1,:));
writetable(T_jasp_mean, sprintf('EMG_Cycles_Moyens_MOYENNES_%s.xlsx', timestamp));

% ========== CALCUL DES RATIOS TRAPÈZE SUPÉRIEUR / SERRATUS ANTÉRIEUR ==========
fprintf('\nCalcul des ratios Trapèze supérieur / Serratus antérieur...\n');

% Identifier les indices des muscles d'intérêt
trapeze_idx = [];
serratus_idx = [];

% Rechercher les muscles (insensible à la casse)
for m = 1:nb_muscles
    muscle_name = muscles_R{m};
    
    % Rechercher Trapèze supérieur
    if contains(muscle_name, 'RTRAPS')
        trapeze_idx = m;
        fprintf('Trapèze supérieur trouvé: %s (index %d)\n', muscles_R{m}, m);
    end
    
    % Rechercher Serratus antérieur
    if contains(muscle_name, 'RSERRA')
        serratus_idx = m;
        fprintf('Serratus antérieur trouvé: %s (index %d)\n', muscles_R{m}, m);
    end
end

% Vérifier que les deux muscles ont été trouvés
if isempty(trapeze_idx) || isempty(serratus_idx)
    fprintf('ERREUR: Impossible de trouver les muscles requis.\n');
    fprintf('Muscles disponibles:\n');
    for m = 1:nb_muscles
        fprintf('  %d: %s\n', m, muscles_R{m});
    end
    fprintf('Veuillez vérifier les noms des muscles.\n');
else
    % Créer le tableau des ratios
    ratio_data = {};
    ratio_data{1, 1} = 'Sujet_ID';
    ratio_data{1, 2} = 'Groupe';
    ratio_data{1, 3} = 'Ratio_Trapeze_Serratus_Moyen';
    
    % Ajouter les ratios pour chaque point temporel du cycle
    for t = 1:normalized_length
        ratio_data{1, 3+t} = sprintf('Ratio_T%d', t);
    end
    
    row_idx = 2;
    sujet_global_id = 1;
    
    % Calculer les ratios pour chaque groupe
    for f = 1:num_files
        group_name = group_labels{f};
        
        % Vérifier que les données existent pour les deux muscles
        if ~isempty(representative_cycles{trapeze_idx, f}) && ~isempty(representative_cycles{serratus_idx, f})
            n_subjects = size(representative_cycles{trapeze_idx, f}, 1);
            
            % Vérifier que le nombre de sujets est cohérent
            if size(representative_cycles{serratus_idx, f}, 1) ~= n_subjects
                fprintf('ATTENTION: Nombre de sujets incohérent entre les muscles pour %s\n', group_name);
                n_subjects = min(size(representative_cycles{trapeze_idx, f}, 1), ...
                               size(representative_cycles{serratus_idx, f}, 1));
            end
            
            for s = 1:n_subjects
                % Extraire les données des deux muscles pour ce sujet
                trapeze_data = representative_cycles{trapeze_idx, f}(s, :);
                serratus_data = representative_cycles{serratus_idx, f}(s, :);
                
                % Calculer les ratios point par point
                ratios_time = trapeze_data ./ serratus_data;
                
                % Remplacer les valeurs infinies ou NaN par NaN
                ratios_time(isinf(ratios_time) | isnan(ratios_time)) = NaN;
                
                % Calculer le ratio moyen sur tout le cycle
                ratio_moyen = mean(ratios_time, 'omitnan');
                
                % Remplir le tableau
                ratio_data{row_idx, 1} = sujet_global_id;
                ratio_data{row_idx, 2} = group_name;
                ratio_data{row_idx, 3} = ratio_moyen;
                
                % Ajouter les ratios pour chaque point temporel
                for t = 1:normalized_length
                    ratio_data{row_idx, 3+t} = ratios_time(t);
                end
                
                row_idx = row_idx + 1;
                sujet_global_id = sujet_global_id + 1;
            end
        else
            fprintf('ATTENTION: Données manquantes pour %s\n', group_name);
        end
    end
    
    % Sauvegarder le tableau des ratios
    ratio_filename = sprintf('Ratios_Trapeze_Serratus_%s.xlsx', timestamp);
    T_ratios = cell2table(ratio_data(2:end, :), 'VariableNames', ratio_data(1, :));
    writetable(T_ratios, ratio_filename, 'Sheet', 'Ratios_Trapeze_Serratus');

    fprintf('Tableau des ratios sauvegardé dans: %s\n', ratio_filename);
    
    % Afficher les statistiques descriptives par groupe
    fprintf('\n=== STATISTIQUES DESCRIPTIVES DES RATIOS ===\n');
    fprintf('%-15s %-10s %-15s %-15s %-15s\n', 'Groupe', 'N', 'Moyenne', 'Écart-type', 'Médiane');
    fprintf(repmat('-', 1, 75));
    fprintf('\n');
    
    for f = 1:num_files
        group_name = group_labels{f};
        
        % Extraire les ratios moyens pour ce groupe
        group_mask = strcmp(T_ratios.Groupe, group_name);
        group_ratios = T_ratios.Ratio_Trapeze_Serratus_Moyen(group_mask);
        
        % Calculer les statistiques
        n_subjects = sum(group_mask);
        mean_ratio = mean(group_ratios, 'omitnan');
        std_ratio = std(group_ratios, 'omitnan');
        median_ratio = median(group_ratios, 'omitnan');
        
        fprintf('%-15s %-10d %-15.3f %-15.3f %-15.3f\n', ...
            group_name, n_subjects, mean_ratio, std_ratio, median_ratio);
    end
    
    fprintf('===============================================\n');
end

fprintf('Calcul des ratios terminé.\n');

    fprintf('\nAnalyse terminée avec succès!\n');
    fprintf('Total de %d muscle(s) analysé(s) avec %d condition(s).\n', nb_muscles, num_files);
    
    if num_files > 1
        fprintf('Figure finale créée avec les cycles représentatifs et les analyses SPM1D.\n');
        fprintf('Les barres horizontales indiquent les périodes significatives pour chaque comparaison.\n');
    end
end