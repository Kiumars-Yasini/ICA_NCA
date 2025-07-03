clc;
close all;
clear;
%%
project='MOBILEM';%, 'MOBILEM','EVERLASTING','MODEL2LIFE', 'cobalt_p'
addpath("DrosteEffect-BrewerMap-3.2.5.0");

%% Main data folder
if strcmp(project,'MOBILEM')
    main_folder = fullfile(pwd,"Mobileem","cycle_mapping_resistance_sorted");
    output_folder = fullfile(pwd,"MOBILEM","ICA_figure");
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
elseif strcmp(project,'EVERLASTING')
    main_folder = fullfile(pwd,"Everlasting","cycle_mapping_resistance_sorted");
    output_folder = fullfile(pwd,"EVERLASTING","ICA_figure");
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
elseif strcmp(project,'MODEL2LIFE')
    main_folder = fullfile(pwd,"Model2Life", "cycle_mapping_resistance_sorted");
    output_folder = fullfile(pwd,"MODEL2LIFE","ICA_figure");
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
elseif strcmp(project,'cobalt_p')
    main_folder = fullfile(pwd,"Cobalt","cycle_mapping_resistance_sorted");
    output_folder = fullfile(pwd,"BMWK_COBALT-P","ICA_figure");
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
end

% Get a list of test condition folders within the main folder
testConditionFolders = dir(main_folder);
testConditionFolders = testConditionFolders([testConditionFolders.isdir]);
testConditionFolders = testConditionFolders(~ismember({testConditionFolders.name}, {'.','..'}));
% -------------------------------------------------------------
%  2) PARFOR: Perform computations in parallel, but DO NOT save
% -------------------------------------------------------------
results_all = struct('cell_struct',{},'testConditionName',{});

for k = 1:length(testConditionFolders)
    % Current test condition
    testConditionName = testConditionFolders(k).name;
    testConditionPath = fullfile(main_folder, testConditionName);

    % Get a list of cell folders
    cellFolders = dir(testConditionPath);
    cellFolders = cellFolders([cellFolders.isdir]);
    cellFolders = cellFolders(~ismember({cellFolders.name}, {'.','..'}));

    % Preallocate a cell to store all results in this test condition
    processed_data = cell(numel(cellFolders), 1);

    for j = 1:length(cellFolders)
        cellFolderPath = fullfile(testConditionPath, cellFolders(j).name);

        % Get a list of .mat files
        matFiles = dir(fullfile(cellFolderPath, '*.mat'));

        for i = 1:length(matFiles)
            matFilePath = fullfile(cellFolderPath, matFiles(i).name);

            % 1) LOAD data into a variable
            data = load(matFilePath); 
            cell_struct = data.cell_struct;  %#ok<NASGU>
            % To do. I want to store the cell_struct in another structure
            results_all(end+1).cell_struct       = cell_struct;
            results_all(end).testConditionName     = testConditionName;

        end
    end
end
          
%%
clc;
for n = 1:length(results_all)
        % Create a local copy of the current cell_struct to work on.
       
        temp_cell_struct = results_all(n).cell_struct;
        thisCond = results_all(n).testConditionName;

        % (Reset computed fields; these will be built up for each step.)
        temp_cell_struct.ICA_x         = {};
        temp_cell_struct.ICA_y         = {};  
        temp_cell_struct.dVA_x         = {};
        temp_cell_struct.dVA_y         = {};
     
        %% Cell chemistry settings
        cell_chemistry_name = 'NCA||GrSi';
        voltage_min = 3.2; 
        voltage_max = 4.8;
        RPT_count  = numel(temp_cell_struct.equivalent_cycle_count);
        equiv   = temp_cell_struct.equivalent_cycle_count(:);
        
        %% Plotting Global style settings
        fontsize = 28;
        linewidth = 2.5;
        interpreter = 'latex';
        setText = @(h, str) set(h, 'String', str, 'FontSize', fontsize, 'Interpreter', interpreter);
        
        % Create and configure figure
        figureHandle = figure('Units', 'inches', 'Position', [0.5, 0.5, 12, 8]);
        set(figureHandle, 'PaperUnits', 'inches', 'PaperSize', [12, 8], ...
                          'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto');
        
        % Make axes on the correct figure
        ax = axes('Parent', figureHandle);
        axesStyle = {'LineWidth', 1.0, ...
                     'XColor', 'k', 'YColor', 'k', ...
                     'TickDir', 'in', 'TickLength', [0.02 0.02], ...
                     'FontSize', fontsize, ...
                     'TickLabelInterpreter', interpreter};
        set(ax, axesStyle{:});
        hold(ax, 'on'); box(ax, 'on');
        
        setText(xlabel(ax, ''), 'Voltage / V');
        setText(ylabel(ax, ''), 'dQ/dV / Ah/V');
        
        cmap = flipud(brewermap(256, 'RdYlBu'));
        colormap(ax, cmap);
        cmin = 0;
        raw_max = temp_cell_struct.equivalent_cycle_count(end);
        if raw_max <=200
            cmax = 200;
            tickVals = cmin:50:cmax;
         elseif (200 <raw_max) && (raw_max <= 500)
            cmax = 500;    
            tickVals = cmin:200:cmax;
         elseif (500 <raw_max) && (raw_max <= 800)
            cmax = 800;
            tickVals = cmin:200:cmax;
         elseif (800 <raw_max) && (raw_max <= 1000)
            cmax = 1000;
            tickVals = cmin:200:cmax;
        elseif (1000 <raw_max) && (raw_max <= 1200)
            cmax = 1200;
            tickVals = cmin:200:cmax;
        elseif (1200 <raw_max) && (raw_max <= 1500)
            cmax = 1500;
            tickVals = cmin:200:cmax;
        elseif (1500 <raw_max) && (raw_max <= 2000)
            cmax = 2000;
            tickVals = cmin:200:cmax;
        else
            cmax = ceil(raw_max / 200) * 200;
            tickVals = cmin:200:cmax;
        end
        
        clim(ax, [cmin cmax]);

        cb = colorbar(ax);
        set(cb, 'Ticks', tickVals, ...
                'TickLabels', tickVals, ...
                'FontSize', fontsize, ...
                'TickLabelInterpreter', interpreter);
        setText(cb.Label, 'Equivalent cycle count');
        xlim(ax, [voltage_min voltage_max]);
        ylim(ax, [0 10 ]);
            % Pre-initialize temporaries so parfor is happy
        Q_smooth   = [];
        dQdV_raw   = [];
        dQdV       = [];

        for l = 1:RPT_count
            Q = temp_cell_struct.AhStep_CHA{1,l}(:);
            U = temp_cell_struct.qOCV_CHA{1,l}(:);
        
        
           %% Actual ICA processing starts here
            
           %% Filter settings
            smoothingMethod = 'butter';
            butterOrder = 4;
            butterFc    = 0.02;
            rloessWin   = 0.5;
            waveletName = 'db8';
            waveletLevel = 3;
            %%
           
           % Filter out non-increasing U
            inc_idx = [true; diff(U) > 0];
            U_filt = U(inc_idx);
            Q_filt = Q(inc_idx);
        
            [Vuniq, idxU] = unique(U_filt);
            Quniq = Q_filt(idxU);
        
            % Smooth
            switch lower(smoothingMethod)
                case 'butter'
                    fs = 1 / mean(diff(Vuniq));
                    [b, a] = butter(butterOrder, butterFc, 'low');
                    Q_smooth = filtfilt(b, a, Quniq);
                case 'rloess'
                    Q_smooth = smoothdata(Quniq, 'rloess', ...
                                          floor(rloessWin * numel(Quniq)));
                case 'wavelet'
                    Q_smooth = wdenoise(Quniq, waveletLevel, ...
                                        'Wavelet', waveletName, ...
                                        'DenoisingMethod', 'SURE');
            end
        
            % Compute dQ/dV
            dQdV_raw = diff(Q_smooth) ./ diff(Vuniq);
            dQdV = filtfilt(ones(1,120)/120, 1, dQdV_raw);
      

            % Get color for current cycle
            xPlot = Vuniq(1:end-1);
            yPlot = dQdV;
            %%
            %Anomaly Detection on Computed ICA Curve ---
        
            %max_dQdV_threshold = 9; 
                                   
           
            %is_outlier = false;
            %if any(yPlot > max_dQdV_threshold)
            %    is_outlier = true;
            %    warning('Cycle %d has an abnormally high dQdV peak (max=%.2f). Flagged as potential outlier.', l, max(yPlot));
            %end

            %if is_outlier
                % You can choose to:
                % 1. Skip plotting this cycle (like before)
                % 2. Store it but mark it (e.g., temp_cell_struct.ICA_y_flagged{end+1} = yPlot)
                % 3. Plot it with a different style (e.g., dashed line, different color)
            %    continue; % This skips plotting and storing for this cycle
            %end
           
            % ICA calculation ends here
            %%
            cv = (equiv(l) - cmin) / (cmax - cmin);
            clr = cmap(max(1, min(256, round(cv*255)+1)), :);
            
            plot(ax, xPlot, yPlot, 'LineWidth', linewidth, 'Color', clr);
                 
            % Append the results for this step to the local (temporary) cell_struct.
            temp_cell_struct.ICA_x{end+1} = xPlot;
            temp_cell_struct.ICA_y{end+1} = yPlot;
        end
        
        % Final plot set up
        cell_name_raw = extractAfter(temp_cell_struct.name, 'M1B_');  % e.g. 'A_1_3'
        cell_name = strrep(cell_name_raw, '_', ' ');             % converts to 'A 1 3'
        title(ax, [cell_chemistry_name, ' cell ID: ', cell_name], 'FontSize', fontsize, 'Interpreter','latex','FontWeight','normal');
        grid(ax, 'on');
        
        % Final rendering and export
        drawnow;
        set(figureHandle, 'Renderer', 'painters');  % Ensure vector output
        % --- Build & ensure final output folder ---
        % TO do for chat gpt here use the respective testConditionFolders.name
        final_output_folder = fullfile(output_folder, thisCond);
        if ~exist(final_output_folder,'dir')
            mkdir(final_output_folder);
        end
            
        % --- Convert any strings to chars ---
        if isstring(cell_name_raw),       cell_name_raw       = char(cell_name_raw);       end
        if isstring(cell_chemistry_name), cell_chemistry_name = char(cell_chemistry_name); end
        if isstring(final_output_folder), final_output_folder = char(final_output_folder); end
        
        % --- Sanitize & build base filename ---
        safe1    = matlab.lang.makeValidName(cell_name_raw);
        safe2    = matlab.lang.makeValidName(cell_chemistry_name);
        file_name = sprintf('%s_%s_ICA', safe1, safe2);   % e.g. 'A_1_8_LFP_ICA'
        
        % --- Full path without extension (char vector) ---
        figure_path = fullfile(final_output_folder, file_name);
        
        % --- Save in multiple formats ---
        % Native .fig
        savefig(figureHandle, [figure_path '.fig']);
        
        % PNG
        saveas(figureHandle, [figure_path '.png']);
        
        % EPS (Level 2 color)
        print(figureHandle, '-depsc2', [figure_path '.eps']);
        
        % SVG
        print(figureHandle, '-dsvg',    [figure_path '.svg']);
        close;
        %% ADD dVA calculation here and 

        %         temp_cell_struct.dVA_x{end+1} = {};
        %         temp_cell_struct.dVA_y{end+1} = {};   

    % Store the updated local structure back into the main array.
    results_all(n).cell_struct = temp_cell_struct;
end

%%

% Get a list of test condition folders within the main folder
testConditionFolders = dir(main_folder);
testConditionFolders = testConditionFolders([testConditionFolders.isdir]);
testConditionFolders = testConditionFolders(~ismember({testConditionFolders.name}, {'.','..'}));

% Initialize an index to iterate through results_all
index = 1;

% Loop through each test condition folder
for k = 1:length(testConditionFolders)
    testConditionName = testConditionFolders(k).name;
    testConditionPath = fullfile(main_folder, testConditionName);

    % Get a list of cell folders within the current test condition folder
    cellFolders = dir(testConditionPath);
    cellFolders = cellFolders([cellFolders.isdir]);
    cellFolders = cellFolders(~ismember({cellFolders.name}, {'.','..'}));

    % Loop through each cell folder
    for j = 1:length(cellFolders)
        cellFolderPath = fullfile(testConditionPath, cellFolders(j).name);

        % Get a list of .mat files in the current cell folder
        matFiles = dir(fullfile(cellFolderPath, '*.mat'));

        % Loop through each .mat file
        for i = 1:length(matFiles)
            matFilePath = fullfile(cellFolderPath, matFiles(i).name);

            % Retrieve the updated cell_struct from results_all
            cell_struct = results_all(index).cell_struct;
           
                        % Extract the folder path from the mat file path
            [destFolder, ~, ~] = fileparts(matFilePath);
            
            % Check if the destination folder exists, if not create it
            if ~exist(destFolder, 'dir')
                mkdir(destFolder);
            end
            % Save the modified cell_struct back to the same .mat file
            save(matFilePath, 'cell_struct');

            % Move to the next element in results_all
            index = index + 1;
        end
    end
end
