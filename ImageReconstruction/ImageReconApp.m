classdef ImageReconApp < matlab.apps.AppBase
    %IMAGERECONAPP First UI version for IMAGE reconstruction workflows.

    properties (Access = public)
        UIFigure matlab.ui.Figure
        Grid matlab.ui.container.GridLayout

        ControlPanel matlab.ui.container.Panel
        ResultPanel matlab.ui.container.Panel
        ParamTabGroup matlab.ui.container.TabGroup
        MainTab matlab.ui.container.Tab
        LseTab matlab.ui.container.Tab
        RenderTab matlab.ui.container.Tab
        AdvancedTab matlab.ui.container.Tab

        InputDirField matlab.ui.control.EditField
        OutputDirField matlab.ui.control.EditField
        LibraryDropDown matlab.ui.control.DropDown
        MethodDropDown matlab.ui.control.DropDown
        ModeDropDown matlab.ui.control.DropDown
        FileModeDropDown matlab.ui.control.DropDown
        BasenameField matlab.ui.control.EditField
        IndicesField matlab.ui.control.EditField
        FileListBox matlab.ui.control.ListBox

        CoarseRadiusField matlab.ui.control.NumericEditField
        TempField matlab.ui.control.NumericEditField
        BaselineField matlab.ui.control.NumericEditField
        TopKField matlab.ui.control.NumericEditField
        BatchSizeField matlab.ui.control.NumericEditField
        SoftmaxCheck matlab.ui.control.CheckBox
        QualityCheck matlab.ui.control.CheckBox
        ResidualCutoffField matlab.ui.control.NumericEditField

        DisplayResolutionField matlab.ui.control.NumericEditField
        GaussianSigmaField matlab.ui.control.NumericEditField
        GaussianRadiusField matlab.ui.control.NumericEditField
        ColorbarModeDropDown matlab.ui.control.DropDown
        ClimHighField matlab.ui.control.NumericEditField
        CdfCheck matlab.ui.control.CheckBox
        UcmCheck matlab.ui.control.CheckBox

        GridPeakRadiusField matlab.ui.control.NumericEditField
        GridSmoothWindowField matlab.ui.control.NumericEditField
        FabbriMaxEventsField matlab.ui.control.EditField
        FabbriStrideField matlab.ui.control.NumericEditField
        FabbriMaxIterField matlab.ui.control.NumericEditField
        AngerPowerField matlab.ui.control.NumericEditField

        RunButton matlab.ui.control.Button
        RerunButton matlab.ui.control.Button
        SelectInputButton matlab.ui.control.Button
        SelectOutputButton matlab.ui.control.Button
        SelectFilesButton matlab.ui.control.Button
        ClearFilesButton matlab.ui.control.Button
        OpenOutputButton matlab.ui.control.Button
        SaveConfigButton matlab.ui.control.Button
        LoadConfigButton matlab.ui.control.Button

        PreviewAxes matlab.ui.control.UIAxes
        SummaryTable matlab.ui.control.Table
        LogTextArea matlab.ui.control.TextArea
    end

    properties (Access = private)
        RepoRoot char
        UiDir char
        CurrentCfg struct
        SelectedFiles cell = {}
        LastResult struct
    end

    methods (Access = private)
        function startupFcn(app)
            app.UiDir = fileparts(mfilename('fullpath'));
            app.RepoRoot = app.UiDir;
            addpath(app.UiDir);
            addpath(fullfile(app.RepoRoot, 'method'));

            app.CurrentCfg = IMAGE_default_config(app.RepoRoot);
            app.loadConfigIntoUi(app.CurrentCfg);
            app.refreshLibraries();
            app.applyModeDefaults();
            app.appendLog('Ready.');
        end

        function refreshLibraries(app)
            libDir = app.CurrentCfg.library_folder;
            files = dir(fullfile(libDir, '*.mat'));
            if isempty(files)
                app.LibraryDropDown.Items = {'No library found'};
                app.LibraryDropDown.Value = 'No library found';
                return;
            end
            names = {files.name};
            app.LibraryDropDown.Items = names;
            if ismember(app.CurrentCfg.lse_library_filename, names)
                app.LibraryDropDown.Value = app.CurrentCfg.lse_library_filename;
            else
                app.LibraryDropDown.Value = names{1};
            end
        end

        function loadConfigIntoUi(app, cfg)
            app.InputDirField.Value = cfg.input_data_folder;
            app.OutputDirField.Value = cfg.output_root;
            app.MethodDropDown.Value = cfg.localization.method;
            app.ModeDropDown.Value = cfg.postprocess_mode;
            app.BasenameField.Value = cfg.file_basename;
            app.IndicesField.Value = app.vectorToText(cfg.file_indices_to_process);

            app.CoarseRadiusField.Value = cfg.localization.coarse_radius;
            app.TempField.Value = cfg.localization.lse_temperature;
            app.BaselineField.Value = cfg.localization.baseline_ratio;
            app.TopKField.Value = cfg.localization.top_k_ratio;
            app.BatchSizeField.Value = cfg.localization.batch_size;
            app.SoftmaxCheck.Value = cfg.localization.enable_softmax_interpolation;
            app.QualityCheck.Value = cfg.localization.enable_quality_filter;
            app.ResidualCutoffField.Value = cfg.localization.residual_percentile_cutoff;

            app.DisplayResolutionField.Value = cfg.display_resolution;
            app.GaussianSigmaField.Value = cfg.render.gaussian_sigma_bins;
            app.GaussianRadiusField.Value = cfg.render.gaussian_radius_bins;
            app.ColorbarModeDropDown.Value = cfg.render.colorbar_mode;
            app.ClimHighField.Value = cfg.render.clim_high_percentile;
            app.CdfCheck.Value = cfg.enable_cdf_correction;
            app.UcmCheck.Value = cfg.enable_ucm_correction;

            app.GridPeakRadiusField.Value = cfg.gridmask.manual_peak_search_radius_mm;
            app.GridSmoothWindowField.Value = cfg.gridmask.profile_smooth_window_bins;
            app.FabbriMaxEventsField.Value = app.scalarToText(cfg.fabbri.max_events_per_file);
            app.FabbriStrideField.Value = cfg.fabbri.event_stride;
            app.FabbriMaxIterField.Value = cfg.fabbri.max_iterations;
            app.AngerPowerField.Value = cfg.anger.projection_power;
        end

        function cfg = collectCfgFromUi(app)
            cfg = app.CurrentCfg;
            cfg.input_data_folder = char(app.InputDirField.Value);
            cfg.output_root = char(app.OutputDirField.Value);
            cfg.localization.method = char(app.MethodDropDown.Value);
            cfg.postprocess_mode = char(app.ModeDropDown.Value);
            cfg.file_basename = char(app.BasenameField.Value);
            cfg.file_indices_to_process = app.parseVector(app.IndicesField.Value);
            cfg.lse_library_filename = char(app.LibraryDropDown.Value);
            cfg.library_file = fullfile(cfg.library_folder, cfg.lse_library_filename);

            cfg.localization.coarse_radius = app.CoarseRadiusField.Value;
            cfg.localization.lse_temperature = app.TempField.Value;
            cfg.localization.baseline_ratio = app.BaselineField.Value;
            cfg.localization.top_k_ratio = app.TopKField.Value;
            cfg.localization.batch_size = app.BatchSizeField.Value;
            cfg.localization.batch_size_64ch = app.BatchSizeField.Value;
            cfg.localization.enable_softmax_interpolation = app.SoftmaxCheck.Value;
            cfg.localization.enable_quality_filter = app.QualityCheck.Value;
            cfg.localization.residual_percentile_cutoff = app.ResidualCutoffField.Value;

            cfg.display_resolution = app.DisplayResolutionField.Value;
            cfg.render.gaussian_sigma_bins = app.GaussianSigmaField.Value;
            cfg.render.gaussian_radius_bins = app.GaussianRadiusField.Value;
            cfg.render.colorbar_mode = char(app.ColorbarModeDropDown.Value);
            cfg.render.clim_high_percentile = app.ClimHighField.Value;
            cfg.enable_cdf_correction = app.CdfCheck.Value;
            cfg.enable_ucm_correction = app.UcmCheck.Value;

            cfg.gridmask.manual_peak_search_radius_mm = app.GridPeakRadiusField.Value;
            cfg.gridmask.profile_smooth_window_bins = app.GridSmoothWindowField.Value;
            cfg.fabbri.max_events_per_file = app.parseScalar(app.FabbriMaxEventsField.Value);
            cfg.fabbri.event_stride = app.FabbriStrideField.Value;
            cfg.fabbri.max_iterations = app.FabbriMaxIterField.Value;
            cfg.anger.projection_power = app.AngerPowerField.Value;

            cfg.selected_files = {};
            if strcmp(app.FileModeDropDown.Value, 'Manual files')
                if isempty(app.SelectedFiles)
                    error('Manual file mode is selected, but no .mat files were chosen.');
                end
                cfg.selected_files = app.SelectedFiles;
            end

            method = strtrim(cfg.localization.method);
            if any(strcmpi(method, {'lse64', 'lse_64ch', 'lse_softmax_64ch'}))
                cfg.localization.batch_size_64ch = cfg.localization.batch_size;
            end
        end

        function runPipeline(app)
            cfg = app.collectCfgFromUi();
            app.CurrentCfg = cfg;
            app.LastResult = struct();
            app.appendLog(sprintf('Running %s / %s ...', cfg.localization.method, cfg.postprocess_mode));

            d = uiprogressdlg(app.UIFigure, 'Title', 'IMAGE reconstruction', ...
                'Message', 'Processing files...', 'Indeterminate', 'on');
            cleanupObj = onCleanup(@() close(d));
            try
                diaryFile = [tempname, '.txt'];
                diary(diaryFile);
                result = IMAGE_run_recon(cfg);
                diary off;
                app.LastResult = result;
                app.showResult(result);
                app.appendDiary(diaryFile);
                app.appendLog('Done.');
            catch e
                diary off;
                app.appendLog(['Error: ', e.message]);
                uialert(app.UIFigure, e.message, 'Run failed');
            end
            clear cleanupObj;
        end

        function showResult(app, result)
            if isfield(result, 'run_summary') && ~isempty(result.run_summary)
                app.SummaryTable.Data = result.run_summary;
            else
                app.SummaryTable.Data = table();
            end

            if isfield(result, 'recon_images') && ~isempty(result.recon_images)
                imagePath = result.recon_images{end};
                if exist(imagePath, 'file')
                    img = imread(imagePath);
                    image(app.PreviewAxes, img);
                    axis(app.PreviewAxes, 'image');
                    app.PreviewAxes.XTick = [];
                    app.PreviewAxes.YTick = [];
                    title(app.PreviewAxes, imagePath, 'Interpreter', 'none');
                end
            end
        end

        function applyModeDefaults(app)
            mode = app.ModeDropDown.Value;
            switch mode
                case 'flood'
                    app.BasenameField.Value = 'calib_events_pflood%04d.mat';
                    app.IndicesField.Value = '0';
                case 'gridmask'
                    app.BasenameField.Value = 'calib_events_pgridmask%04d.mat';
                    app.IndicesField.Value = '1010';
                case 'slit'
                    app.BasenameField.Value = 'calib_events_pslity%04d.mat';
                    app.IndicesField.Value = '0';
                case 'single_hole'
                    app.BasenameField.Value = 'calib_events_p%04d.mat';
                    app.IndicesField.Value = '801:809';
                case 'single_hole_scan'
                    app.BasenameField.Value = 'calib_events_p%04d.mat';
                    app.IndicesField.Value = '801:809';
            end
        end

        function MethodChanged(app, ~)
            method = app.MethodDropDown.Value;
            if strcmp(method, 'lse_softmax_64ch') && app.BatchSizeField.Value > 1000
                app.BatchSizeField.Value = 1000;
                app.appendLog('64-channel LSE selected: batch size set to 1000.');
            end
            if strcmp(method, 'anger_standard')
                app.CurrentCfg.anger.position_weight_mode = 'standard';
            elseif strcmp(method, 'anger_RTP')
                app.CurrentCfg.anger.position_weight_mode = 'rtp';
            end
        end

        function FileModeChanged(app, ~)
            isManual = strcmp(app.FileModeDropDown.Value, 'Manual files');
            app.SelectFilesButton.Enable = app.boolToOnOff(isManual);
            app.ClearFilesButton.Enable = app.boolToOnOff(isManual);
            app.BasenameField.Enable = app.boolToOnOff(~isManual);
            app.IndicesField.Enable = app.boolToOnOff(~isManual);
        end

        function SelectInputDir(app, ~)
            pathName = uigetdir(app.InputDirField.Value, 'Select input data folder');
            if isequal(pathName, 0)
                return;
            end
            app.InputDirField.Value = pathName;
        end

        function SelectOutputDir(app, ~)
            pathName = uigetdir(app.OutputDirField.Value, 'Select output folder');
            if isequal(pathName, 0)
                return;
            end
            app.OutputDirField.Value = pathName;
        end

        function SelectFiles(app, ~)
            [files, pathName] = uigetfile('*.mat', 'Select planeset .mat files', ...
                app.InputDirField.Value, 'MultiSelect', 'on');
            if isequal(files, 0)
                return;
            end
            if ischar(files)
                files = {files};
            end
            app.SelectedFiles = cellfun(@(name) fullfile(pathName, name), files, 'UniformOutput', false);
            app.FileListBox.Items = app.SelectedFiles;
            app.appendLog(sprintf('Selected %d file(s).', numel(app.SelectedFiles)));
        end

        function ClearFiles(app, ~)
            app.SelectedFiles = {};
            app.FileListBox.Items = {};
        end

        function OpenOutput(app, ~)
            outputRoot = app.OutputDirField.Value;
            if ~exist(outputRoot, 'dir')
                uialert(app.UIFigure, 'Output folder does not exist yet.', 'Open output');
                return;
            end
            winopen(outputRoot);
        end

        function SaveConfig(app, ~)
            cfg = app.collectCfgFromUi();
            [file, pathName] = uiputfile('*.mat', 'Save IMAGE UI config', ...
                fullfile(app.UiDir, 'image_ui_config.mat'));
            if isequal(file, 0)
                return;
            end
            save(fullfile(pathName, file), 'cfg');
            app.appendLog(['Saved config: ', fullfile(pathName, file)]);
        end

        function LoadConfig(app, ~)
            [file, pathName] = uigetfile('*.mat', 'Load IMAGE UI config', app.UiDir);
            if isequal(file, 0)
                return;
            end
            data = load(fullfile(pathName, file), 'cfg');
            if ~isfield(data, 'cfg')
                uialert(app.UIFigure, 'Selected file does not contain cfg.', 'Load config');
                return;
            end
            app.CurrentCfg = data.cfg;
            app.loadConfigIntoUi(app.CurrentCfg);
            app.refreshLibraries();
            app.appendLog(['Loaded config: ', fullfile(pathName, file)]);
        end

        function appendDiary(app, diaryFile)
            if exist(diaryFile, 'file')
                txt = fileread(diaryFile);
                if ~isempty(strtrim(txt))
                    app.appendLog(txt);
                end
                delete(diaryFile);
            end
        end

        function appendLog(app, msg)
            if isempty(msg)
                return;
            end
            msg = string(msg);
            parts = splitlines(msg);
            parts(parts == "") = [];
            if isempty(parts)
                return;
            end
            stamp = string(datetime('now', 'Format', 'HH:mm:ss'));
            lines = cellstr(app.LogTextArea.Value);
            for ii = 1:numel(parts)
                lines{end+1} = char(stamp + "  " + parts(ii)); %#ok<AGROW>
            end
            if numel(lines) > 300
                lines = lines(end-299:end);
            end
            app.LogTextArea.Value = lines;
            drawnow limitrate;
        end

        function indices = parseVector(~, textValue)
            textValue = strtrim(char(textValue));
            if isempty(textValue)
                indices = [];
                return;
            end
            indices = str2num(textValue); %#ok<ST2NM>
            if isempty(indices)
                error('Invalid index expression: %s', textValue);
            end
            indices = reshape(indices, 1, []);
        end

        function value = parseScalar(~, textValue)
            textValue = strtrim(char(textValue));
            if strcmpi(textValue, 'inf')
                value = Inf;
                return;
            end
            value = str2double(textValue);
            if ~isfinite(value)
                error('Invalid scalar value: %s', textValue);
            end
        end

        function textValue = vectorToText(~, values)
            if isempty(values)
                textValue = '';
            elseif numel(values) > 2 && all(diff(values) == diff(values(1:2)))
                step = diff(values(1:2));
                if step == 1
                    textValue = sprintf('%g:%g', values(1), values(end));
                else
                    textValue = sprintf('%g:%g:%g', values(1), step, values(end));
                end
            else
                textValue = mat2str(values);
            end
        end

        function textValue = scalarToText(~, value)
            if isinf(value)
                textValue = 'Inf';
            else
                textValue = num2str(value);
            end
        end

        function value = boolToOnOff(~, flag)
            if flag
                value = 'on';
            else
                value = 'off';
            end
        end
    end

    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Name', 'IMAGE Reconstruction UI');
            app.UIFigure.Position = [80 80 1240 760];

            app.Grid = uigridlayout(app.UIFigure, [1 2]);
            app.Grid.ColumnWidth = {430, '1x'};
            app.Grid.RowHeight = {'1x'};

            app.ControlPanel = uipanel(app.Grid, 'Title', 'Controls');
            app.ControlPanel.Layout.Row = 1;
            app.ControlPanel.Layout.Column = 1;

            controlGrid = uigridlayout(app.ControlPanel, [5 1]);
            controlGrid.RowHeight = {98, 145, '1x', 40, 40};
            controlGrid.Padding = [10 8 10 8];

            pathGrid = uigridlayout(controlGrid, [3 3]);
            pathGrid.ColumnWidth = {70, '1x', 78};
            pathGrid.RowHeight = {26, 26, 26};
            uilabel(pathGrid, 'Text', 'Input');
            app.InputDirField = uieditfield(pathGrid, 'text');
            app.InputDirField.Layout.Row = 1;
            app.InputDirField.Layout.Column = 2;
            app.SelectInputButton = uibutton(pathGrid, 'Text', 'Browse', 'ButtonPushedFcn', @(src,event) app.SelectInputDir(event));
            app.SelectInputButton.Layout.Row = 1;
            app.SelectInputButton.Layout.Column = 3;
            uilabel(pathGrid, 'Text', 'Output');
            app.OutputDirField = uieditfield(pathGrid, 'text');
            app.OutputDirField.Layout.Row = 2;
            app.OutputDirField.Layout.Column = 2;
            app.SelectOutputButton = uibutton(pathGrid, 'Text', 'Browse', 'ButtonPushedFcn', @(src,event) app.SelectOutputDir(event));
            app.SelectOutputButton.Layout.Row = 2;
            app.SelectOutputButton.Layout.Column = 3;
            uilabel(pathGrid, 'Text', 'Library');
            app.LibraryDropDown = uidropdown(pathGrid);
            app.LibraryDropDown.Layout.Row = 3;
            app.LibraryDropDown.Layout.Column = [2 3];

            topGrid = uigridlayout(controlGrid, [4 2]);
            topGrid.ColumnWidth = {120, '1x'};
            topGrid.RowHeight = {26, 26, 26, 26};
            uilabel(topGrid, 'Text', 'Method');
            app.MethodDropDown = uidropdown(topGrid, 'Items', ...
                {'lse_softmax_64ch', 'lse_softmax', 'fabbri_analytical_lsf', 'anger_standard', 'anger_RTP'}, ...
                'ValueChangedFcn', @(src,event) app.MethodChanged(event));
            uilabel(topGrid, 'Text', 'Postprocess');
            app.ModeDropDown = uidropdown(topGrid, 'Items', ...
                {'flood', 'gridmask', 'slit', 'single_hole', 'single_hole_scan'}, ...
                'ValueChangedFcn', @(src,event) app.applyModeDefaults());
            uilabel(topGrid, 'Text', 'File mode');
            app.FileModeDropDown = uidropdown(topGrid, 'Items', {'Pattern + indices', 'Manual files'}, ...
                'ValueChangedFcn', @(src,event) app.FileModeChanged(event));
            uilabel(topGrid, 'Text', 'Basename');
            app.BasenameField = uieditfield(topGrid, 'text');

            app.ParamTabGroup = uitabgroup(controlGrid);
            app.MainTab = uitab(app.ParamTabGroup, 'Title', 'Files');
            app.LseTab = uitab(app.ParamTabGroup, 'Title', 'Method');
            app.RenderTab = uitab(app.ParamTabGroup, 'Title', 'Render');
            app.AdvancedTab = uitab(app.ParamTabGroup, 'Title', 'Advanced');

            fileGrid = uigridlayout(app.MainTab, [5 2]);
            fileGrid.ColumnWidth = {120, '1x'};
            fileGrid.RowHeight = {26, 32, 32, '1x', 26};
            uilabel(fileGrid, 'Text', 'Indices');
            app.IndicesField = uieditfield(fileGrid, 'text');
            app.IndicesField.Layout.Column = 2;
            app.SelectFilesButton = uibutton(fileGrid, 'Text', 'Select .mat files', ...
                'Enable', 'off', 'ButtonPushedFcn', @(src,event) app.SelectFiles(event));
            app.SelectFilesButton.Layout.Row = 2;
            app.SelectFilesButton.Layout.Column = [1 2];
            app.ClearFilesButton = uibutton(fileGrid, 'Text', 'Clear selected files', ...
                'Enable', 'off', 'ButtonPushedFcn', @(src,event) app.ClearFiles(event));
            app.ClearFilesButton.Layout.Row = 3;
            app.ClearFilesButton.Layout.Column = [1 2];
            app.FileListBox = uilistbox(fileGrid);
            app.FileListBox.Layout.Row = 4;
            app.FileListBox.Layout.Column = [1 2];
            uilabel(fileGrid, 'Text', 'Manual mode ignores basename/indices.');

            methodGrid = uigridlayout(app.LseTab, [8 2]);
            methodGrid.ColumnWidth = {175, '1x'};
            methodGrid.RowHeight = repmat({27}, 1, 8);
            uilabel(methodGrid, 'Text', 'Coarse radius');
            app.CoarseRadiusField = uieditfield(methodGrid, 'numeric');
            uilabel(methodGrid, 'Text', 'LSE temperature');
            app.TempField = uieditfield(methodGrid, 'numeric');
            uilabel(methodGrid, 'Text', 'Baseline ratio');
            app.BaselineField = uieditfield(methodGrid, 'numeric');
            uilabel(methodGrid, 'Text', 'Top-k ratio');
            app.TopKField = uieditfield(methodGrid, 'numeric');
            uilabel(methodGrid, 'Text', 'Batch size');
            app.BatchSizeField = uieditfield(methodGrid, 'numeric');
            app.SoftmaxCheck = uicheckbox(methodGrid, 'Text', 'Softmax interpolation');
            app.SoftmaxCheck.Layout.Column = [1 2];
            app.QualityCheck = uicheckbox(methodGrid, 'Text', 'Quality filter');
            app.QualityCheck.Layout.Column = [1 2];
            uilabel(methodGrid, 'Text', 'Residual cutoff pct');
            app.ResidualCutoffField = uieditfield(methodGrid, 'numeric');

            renderGrid = uigridlayout(app.RenderTab, [7 2]);
            renderGrid.ColumnWidth = {175, '1x'};
            renderGrid.RowHeight = repmat({27}, 1, 7);
            uilabel(renderGrid, 'Text', 'Display resolution');
            app.DisplayResolutionField = uieditfield(renderGrid, 'numeric');
            uilabel(renderGrid, 'Text', 'Gaussian sigma bins');
            app.GaussianSigmaField = uieditfield(renderGrid, 'numeric');
            uilabel(renderGrid, 'Text', 'Gaussian radius bins');
            app.GaussianRadiusField = uieditfield(renderGrid, 'numeric');
            uilabel(renderGrid, 'Text', 'Colorbar mode');
            app.ColorbarModeDropDown = uidropdown(renderGrid, 'Items', {'counts_per_mm2', 'counts'});
            uilabel(renderGrid, 'Text', 'High percentile');
            app.ClimHighField = uieditfield(renderGrid, 'numeric');
            app.CdfCheck = uicheckbox(renderGrid, 'Text', 'CDF correction');
            app.CdfCheck.Layout.Column = [1 2];
            app.UcmCheck = uicheckbox(renderGrid, 'Text', 'UCM correction');
            app.UcmCheck.Layout.Column = [1 2];

            advGrid = uigridlayout(app.AdvancedTab, [6 2]);
            advGrid.ColumnWidth = {175, '1x'};
            advGrid.RowHeight = repmat({27}, 1, 6);
            uilabel(advGrid, 'Text', 'Grid peak radius mm');
            app.GridPeakRadiusField = uieditfield(advGrid, 'numeric');
            uilabel(advGrid, 'Text', 'Grid smooth window');
            app.GridSmoothWindowField = uieditfield(advGrid, 'numeric');
            uilabel(advGrid, 'Text', 'Fabbri max events');
            app.FabbriMaxEventsField = uieditfield(advGrid, 'text');
            uilabel(advGrid, 'Text', 'Fabbri event stride');
            app.FabbriStrideField = uieditfield(advGrid, 'numeric');
            uilabel(advGrid, 'Text', 'Fabbri max iterations');
            app.FabbriMaxIterField = uieditfield(advGrid, 'numeric');
            uilabel(advGrid, 'Text', 'Anger projection power');
            app.AngerPowerField = uieditfield(advGrid, 'numeric');

            runGrid = uigridlayout(controlGrid, [1 2]);
            app.RunButton = uibutton(runGrid, 'Text', 'Run', 'ButtonPushedFcn', @(src,event) app.runPipeline());
            app.RerunButton = uibutton(runGrid, 'Text', 'Rerun current', 'ButtonPushedFcn', @(src,event) app.runPipeline());

            utilGrid = uigridlayout(controlGrid, [1 3]);
            app.OpenOutputButton = uibutton(utilGrid, 'Text', 'Open output', 'ButtonPushedFcn', @(src,event) app.OpenOutput(event));
            app.SaveConfigButton = uibutton(utilGrid, 'Text', 'Save cfg', 'ButtonPushedFcn', @(src,event) app.SaveConfig(event));
            app.LoadConfigButton = uibutton(utilGrid, 'Text', 'Load cfg', 'ButtonPushedFcn', @(src,event) app.LoadConfig(event));

            app.ResultPanel = uipanel(app.Grid, 'Title', 'Results');
            app.ResultPanel.Layout.Row = 1;
            app.ResultPanel.Layout.Column = 2;
            resultGrid = uigridlayout(app.ResultPanel, [3 1]);
            resultGrid.RowHeight = {'1.4x', '0.8x', 150};

            app.PreviewAxes = uiaxes(resultGrid);
            title(app.PreviewAxes, 'Preview');
            app.SummaryTable = uitable(resultGrid);
            app.LogTextArea = uitextarea(resultGrid, 'Editable', 'off');
        end
    end

    methods (Access = public)
        function app = ImageReconApp
            createComponents(app);
            registerApp(app, app.UIFigure);
            runStartupFcn(app, @startupFcn);
            if nargout == 0
                clear app;
            end
        end

        function delete(app)
            delete(app.UIFigure);
        end
    end
end
