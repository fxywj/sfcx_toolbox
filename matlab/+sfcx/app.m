function fig = app(lang)
%APP SFCX graphical front-end: import, variability analysis, coupling values.
%
%   sfcx.app          % Chinese UI
%   sfcx.app('en')    % English UI
%
%   Tab 1  数据导入    single-subject or multi-subject batch import
%   Tab 2  变异分析    dispersion axes, layers, epicentres, sweep  (n >= 2)
%   Tab 3  耦合值计算  per-subject coupling magnitudes + export  (n >= 1)
%   Tab 4  日志        run log and variability-result export
if nargin < 1 || isempty(lang)
    lang = getpref('sfcx', 'lang', 'zh');
end

S = struct();
S.lang = lang;
S.mode = 'single';
S.scSingle = []; S.fcSingle = []; S.boldSingle = [];
S.scBatch = []; S.fcBatch = []; S.boldBatch = [];
S.networks = []; S.hierarchy = []; S.coords = [];
S.behaviour = []; S.covariates = [];
S.coh = []; S.bundle = []; S.vals = [];
S.demoMode = false;
S.h = struct();
S.status = struct();

fig = uifigure('Name', T('appName'), 'Position', [40 40 1320 850]);
fig.Color = [0.96 0.97 0.98];

root = uigridlayout(fig, [1 1]);
root.Padding = [10 10 10 10];
tabs = uitabgroup(root);
S.h.tabs = tabs;
S.h.tData = uitab(tabs);   S.h.tVar = uitab(tabs);
S.h.tMag = uitab(tabs);    S.h.tLog = uitab(tabs);

buildDataTab();
buildVarTab();
buildMagTab();
buildLogTab();
applyLang();
S.progress.Text = T('ready');
logLine(sfcx.describe());
logLine(T('logHead'));

% ===================================================================
%  TAB 1  DATA
% ===================================================================
    function buildDataTab()
        g = uigridlayout(S.h.tData, [2 3]);
        g.RowHeight = {'1x', 26};
        g.ColumnWidth = {300, 430, '1x'};
        g.ColumnSpacing = 12;
        g.Padding = [6 6 6 6];

        % ---- left: mode + build -----------------------------------
        S.h.pMode = uipanel(g, 'BackgroundColor', [1 1 1]);
        S.h.pMode.Layout.Row = 1;
        S.h.pMode.Layout.Column = 1;
        gm = uigridlayout(S.h.pMode, [5 1]);
        gm.RowHeight = {64, 58, 30, 8, 52};
        gm.Padding = [8 8 8 8];
        S.h.modeGroup = uibuttongroup(gm, 'BorderType', 'none', ...
            'SelectionChangedFcn', @onModeChange);
        S.h.modeGroup.Layout.Row = 1;
        S.h.rb1 = uiradiobutton(S.h.modeGroup, 'Text', T('modeSingle'), ...
                                'Value', true, 'FontWeight', 'bold');
        S.h.rb1.Position = [12 40 260 20];
        S.h.rb2 = uiradiobutton(S.h.modeGroup, 'Text', T('modeMulti'), ...
                                'FontWeight', 'bold');
        S.h.rb2.Position = [12 14 260 20];
        S.h.langBtn = uibutton(gm, 'push', 'ButtonPushedFcn', @onLang);
        S.h.langBtn.Layout.Row = 3;
        S.h.modeHint = uilabel(gm, 'FontColor', [0.42 0.49 0.55], ...
                               'WordWrap', 'on');
        S.h.modeHint.Layout.Row = 2;
        S.h.modeHint.Layout.Column = 1;
        S.h.buildBtn = uibutton(gm, 'push', 'FontWeight', 'bold', ...
            'BackgroundColor', [0.16 0.42 0.36], 'FontColor', [1 1 1], ...
            'ButtonPushedFcn', @onBuild);
        S.h.buildBtn.Layout.Row = 4;

        % ---- middle: import slots (single & batch share the panel) --
        S.h.pSlot = uipanel(g, 'Title', T('secSingle'), 'BackgroundColor', [1 1 1]);
        S.h.pSlot.Layout.Row = 1;
        S.h.pSlot.Layout.Column = 2;
        gsl = uigridlayout(S.h.pSlot, [8 2]);
        gsl.ColumnWidth = {110, '1x'};
        gsl.RowHeight = {34, 34, 34, 34, 34, 34, 34, 34};
        gsl.RowSpacing = 7;
        gsl.Padding = [8 8 8 8];
        S.h.slotGrid = gsl;
        addSlot(gsl, 1, 'sc', @onLoadSCSingle, 'scLab');
        addSlot(gsl, 2, 'fc', @onLoadFCSingle, 'fcLab');
        addSlot(gsl, 3, 'bold', @onLoadBoldSingle, 'boldLab');
        S.h.demoBtn = uibutton(gsl, 'push', 'ButtonPushedFcn', @onDemo);
        S.h.demoBtn.Layout.Row = 4;
        S.h.demoBtn.Layout.Column = [1 2];
        addSlot(gsl, 5, 'scBatch', @onLoadSCBatch, 'scBatchLab');
        addSlot(gsl, 6, 'fcBatch', @onLoadFCBatch, 'fcBatchLab');
        addSlot(gsl, 7, 'boldBatch', @onLoadBoldBatch, 'boldBatchLab');
        S.h.annoBtn = uibutton(gsl, 'push', 'ButtonPushedFcn', @onLoadAnno);
        S.h.annoBtn.Layout.Row = 8;
        S.h.annoBtn.Layout.Column = [1 2];

        % ---- right: annotation slots -------------------------------
        S.h.pBuild = uipanel(g, 'Title', T('secBuild'), 'BackgroundColor', [1 1 1]);
        S.h.pBuild.Layout.Row = 1;
        S.h.pBuild.Layout.Column = 3;
        gp = uigridlayout(S.h.pBuild, [5 2]);
        gp.ColumnWidth = {96, '1x'};
        gp.RowHeight = {34, 34, 34, 34, 34};
        gp.RowSpacing = 7;
        gp.Padding = [8 8 8 8];
        addAnnoRow(gp, 1, 'networks', @onLoadNetworks, 'netLab');
        addAnnoRow(gp, 2, 'hierarchy', @onLoadHierarchy, 'hierLab');
        addAnnoRow(gp, 3, 'coords', @onLoadCoords, 'coordLab');
        addAnnoRow(gp, 4, 'behaviour', @onLoadBehaviour, 'behLab');
        addAnnoRow(gp, 5, 'covariates', @onLoadCovariates, 'covLab');

        S.progress = uilabel(g, 'FontColor', [0.42 0.49 0.55]);
        S.progress.Layout.Row = 2;
        S.progress.Layout.Column = [1 3];
        updateModePanels();
    end

    function addSlot(g, row, key, cb, nameKey)
        nameLab = uilabel(g, 'Text', T(nameKey), 'HorizontalAlignment', 'left');
        nameLab.Layout.Row = row;
        nameLab.Layout.Column = 1;
        S.h.(['slotLab' key]) = nameLab;
        rowG = uigridlayout(g);
        rowG.Layout.Row = row;
        rowG.Layout.Column = 2;
        rowG.ColumnWidth = {96, '1x'};
        rowG.RowHeight = {30};
        rowG.Padding = [0 0 0 0];
        btn = uibutton(rowG, 'push', 'Text', T('load'), 'ButtonPushedFcn', cb);
        btn.Layout.Column = 1;
        lab = uilabel(rowG, 'Text', T('empty'), 'FontColor', [0.42 0.49 0.55]);
        lab.Layout.Column = 2;
        S.h.(['btn' key]) = btn;
        S.status.(key) = lab;
    end

    function addAnnoRow(g, row, key, cb, nameKey)
        nameLab = uilabel(g, 'Text', T(nameKey), 'HorizontalAlignment', 'left');
        nameLab.Layout.Row = row;
        nameLab.Layout.Column = 1;
        S.h.(['slotLab' key]) = nameLab;
        rowG = uigridlayout(g);
        rowG.Layout.Row = row;
        rowG.Layout.Column = 2;
        rowG.ColumnWidth = {80, '1x'};
        rowG.RowHeight = {28};
        rowG.Padding = [0 0 0 0];
        btn = uibutton(rowG, 'push', 'Text', T('load'), 'ButtonPushedFcn', cb);
        btn.Layout.Column = 1;
        lab = uilabel(rowG, 'Text', T('empty'), 'FontColor', [0.42 0.49 0.55]);
        lab.Layout.Column = 2;
        S.h.(['btnAnno_' key]) = btn;
        S.status.(key) = lab;
    end

    function onLang(~, ~)
        if strcmp(S.lang, 'zh')
            S.lang = 'en';
        else
            S.lang = 'zh';
        end
        setpref('sfcx', 'lang', S.lang);
        applyLang();
        S.progress.Text = T('ready');
    end

    function onModeChange(src, ~)
        if src.SelectedObject == S.h.rb1
            S.mode = 'single';
        else
            S.mode = 'multi';
        end
        updateModePanels();
    end

    function updateModePanels()
        gl = S.h.slotGrid;
        if strcmp(S.mode, 'single')
            S.h.pSlot.Title = T('secSingle');
            S.h.modeHint.Text = T('modeSingleHint');
            gl.RowHeight = {34, 34, 34, 34, 1, 1, 1, 1};
        else
            S.h.pSlot.Title = T('secBatch');
            S.h.modeHint.Text = T('modeMultiHint');
            gl.RowHeight = {1, 1, 1, 1, 34, 34, 34, 34};
        end
        vis = strcmp(S.mode, 'single');
        S.h.btnsc.Visible = visOn(vis);   S.h.btnfc.Visible = visOn(vis);
        S.h.btnbold.Visible = visOn(vis); S.h.demoBtn.Visible = visOn(vis);
        S.h.btnscBatch.Visible = visOn(~vis);
        S.h.btnfcBatch.Visible = visOn(~vis);
        S.h.btnboldBatch.Visible = visOn(~vis);
        S.h.annoBtn.Visible = 'on';
        S.status.sc.Visible = visOn(vis);       S.status.fc.Visible = visOn(vis);
        S.status.bold.Visible = visOn(vis);     S.status.scBatch.Visible = visOn(~vis);
        S.status.fcBatch.Visible = visOn(~vis); S.status.boldBatch.Visible = visOn(~vis);
        S.h.slotLabSc.Visible = visOn(vis);     S.h.slotLabFc.Visible = visOn(vis);
        S.h.slotLabBold.Visible = visOn(vis);
        S.h.slotLabScBatch.Visible = visOn(~vis);
        S.h.slotLabFcBatch.Visible = visOn(~vis);
        S.h.slotLabBoldBatch.Visible = visOn(~vis);
    end

    function s = visOn(tf)
        if tf
            s = 'on';
        else
            s = 'off';
        end
    end

% ------------------------------- data callbacks --------------------
    function onLoadSCSingle(~, ~), S.scSingle = pickAny(T('pickSC'), 'sc'); end
    function onLoadFCSingle(~, ~), S.fcSingle = pickAny(T('pickFC'), 'fc'); end
    function onLoadBoldSingle(~, ~), S.boldSingle = pickAny(T('pickBold'), 'bold'); end

    function onLoadSCBatch(~, ~), S.scBatch = pickBatch('sc'); end
    function onLoadFCBatch(~, ~), S.fcBatch = pickBatch('fc'); end
    function onLoadBoldBatch(~, ~), S.boldBatch = pickBatch('bold'); end

    function onLoadAnno(~, ~)
        s = loadAny(T('pickAnno'));
        if isempty(s)
            return
        end
        filled = grabAnno(s, {'networks', 'hierarchy', 'coords', 'behaviour', 'covariates'});
        if isempty(filled)
            uialert(fig, T('annoNoneLong'), T('annoLab'));
        else
            S.progress.Text = sprintf(T('smartFilled'), strjoin(filled, ', '));
        end
    end

    function onLoadNetworks(~, ~), pickAnnoOne('networks'); end
    function onLoadHierarchy(~, ~), pickAnnoOne('hierarchy'); end
    function onLoadCoords(~, ~), pickAnnoOne('coords'); end
    function onLoadBehaviour(~, ~), pickAnnoOne('behaviour'); end
    function onLoadCovariates(~, ~), pickAnnoOne('covariates'); end

    function pickAnnoOne(key)
        s = loadAny(T('pickAnno'));
        if isempty(s)
            return
        end
        if isfield(s, key)
            S.(key) = s.(key);
            S.status.(key).Text = sizeText(S.(key));
        else
            uialert(fig, sprintf('%s: %s', key, strjoin(fieldnames(s)', ', ')), T('loadTitle'));
        end
    end

    function filled = grabAnno(s, keys)
        filled = {};
        for i = 1:numel(keys)
            k = keys{i};
            if isfield(s, k) && ~isempty(S.(k))
                continue
            end
            if isfield(s, k)
                S.(k) = s.(k);
                S.status.(k).Text = sizeText(S.(k));
                filled{end + 1} = k; %#ok<AGROW>
            end
        end
    end

    function onDemo(~, ~)
        S.coh = sfcx.demo(160, 100, 0);
        S.demoMode = true;
        setSlotText('sc', 'demo 160×100');
        setSlotText('fc', 'demo 160×100');
        setSlotText('bold', T('empty'));
        setSlotText('scBatch', 'demo 160×100');
        setSlotText('fcBatch', 'demo 160×100');
        setSlotText('boldBatch', T('empty'));
        S.progress.Text = T('demoReady');
        logLine('demo: 160 subjects x 100 regions');
    end

    function onBuild(~, ~)
        try
            if S.demoMode && ~isempty(S.coh)
                showCohort(S.coh);
                return
            end
            cfg = struct();
            if strcmp(S.mode, 'single')
                cfg.sc = S.scSingle;
                cfg.fc = S.fcSingle;
                cfg.bold = S.boldSingle;
            else
                cfg.sc = S.scBatch;
                cfg.fc = S.fcBatch;
                cfg.bold = S.boldBatch;
            end
            if isempty(cfg.sc) && isempty(cfg.fc) && isempty(cfg.bold)
                uialert(fig, T('needData'), T('secBuild'));
                return
            end
            cfg.networks = S.networks;
            cfg.hierarchy = S.hierarchy;
            cfg.coords = S.coords;
            cfg.behaviour = S.behaviour;
            cfg.covariates = S.covariates;
            [coh, rep] = sfcx.buildCohort(cfg);
            S.coh = coh;
            S.bundle = [];
            S.vals = [];
            S.demoMode = false;
            showCohort(coh);
            for i = 1:numel(rep.warnings)
                logLine(['! ' rep.warnings{i}]);
            end
        catch err
            uialert(fig, err.message, T('buildFail'));
        end
    end

    function showCohort(coh)
        S.progress.Text = sprintf(T('builtFmt'), coh.n, coh.nNodes, coh.nEdges);
        if coh.n < 2
            logLine(sprintf('!! %s', sprintf(T('oneSubjectWarn'), coh.n)));
        elseif coh.n < 8
            logLine(sprintf('!! %s', sprintf(T('smallWarn'), coh.n)));
        end
        logLine(sprintf('cohort: %d subjects x %d nodes x %d edges', ...
                        coh.n, coh.nNodes, coh.nEdges));
    end

% ------------------------------- loaders ---------------------------
    function v = pickAny(prompt, kind)
        v = [];
        [file, path] = uigetfile({'*.mat;*.csv;*.txt;*.tsv', 'data'}, prompt);
        if isequal(file, 0)
            return
        end
        full = fullfile(path, file);
        try
            v = readOne(full, kind);
            setSlotText(kind, sprintf('%s  %s', file, sizeText(v)));
        catch err
            uialert(fig, err.message, T('loadTitle'));
        end
    end

    function v = pickBatch(kind)
        v = [];
        choice = uiconfirm(fig, T('batchChoice'), T('batchTitle'), ...
                           'Options', {T('batchFolder'), T('batchMat'), T('cancel')});
        if strcmp(choice, T('cancel'))
            return
        end
        try
            if strcmp(choice, T('batchFolder'))
                d = uigetdir(pwd, T('pickBatchDir'));
                if isequal(d, 0)
                    return
                end
                [v, ids] = sfcx.batchStack(d, kind);
                setSlotText(kind, sprintf('%d files  %s', numel(ids), sizeText(v)));
            else
                [file, path] = uigetfile({'*.mat', 'MAT'}, T('pickBatchMat'));
                if isequal(file, 0)
                    return
                end
                v = matTo3D(fullfile(path, file));
                setSlotText(kind, sprintf('%s  %s', file, sizeText(v)));
            end
        catch err
            uialert(fig, err.message, T('loadTitle'));
        end
    end

    function setSlotText(key, txt)
        if isfield(S.status, key)
            S.status.(key).Text = txt;
        end
    end

% ===================================================================
%  TAB 2  VARIABILITY
% ===================================================================
    function buildVarTab()
        g = uigridlayout(S.h.tVar, [1 2]);
        g.ColumnWidth = {320, '1x'};
        g.ColumnSpacing = 12;
        g.Padding = [6 6 6 6];

        S.h.pVarOpts = uipanel(g, 'Title', T('secVarOpts'), 'BackgroundColor', [1 1 1]);
        S.h.pVarOpts.Layout.Column = 1;
        gv = uigridlayout(S.h.pVarOpts, [7 2]);
        gv.ColumnWidth = {86, '1x'};
        gv.RowHeight = {18, 170, 28, 28, 28, 12, 40};
        gv.RowSpacing = 6;
        gv.Padding = [8 8 8 8];

        S.h.varMetricHint = uilabel(gv, 'FontColor', [0.42 0.49 0.55]);
        S.h.varMetricHint.Layout.Row = 1;
        S.h.varMetricHint.Layout.Column = [1 2];
        S.varMetrics = uilistbox(gv, 'Multiselect', 'on');
        S.varMetrics.Layout.Row = 2;
        S.varMetrics.Layout.Column = [1 2];
        S.h.labEstV = uilabel(gv, 'Text', T('estimator'), 'HorizontalAlignment', 'right');
        S.h.labEstV.Layout.Row = 3;
        S.estDrop = uidropdown(gv, 'Items', {'all', 'regional', 'node', 'edge'}, ...
                               'ItemsData', {'all', 'regional', 'node', 'edge'}, ...
                               'Value', 'all');
        S.estDrop.Layout.Row = 3;
        S.estDrop.Layout.Column = 2;
        S.h.labNullV = uilabel(gv, 'Text', T('nullModel'), 'HorizontalAlignment', 'right');
        S.h.labNullV.Layout.Row = 4;
        S.nullDrop = uidropdown(gv, 'Items', {'auto', 'spin', 'msr', 'permutation'}, ...
                                'ItemsData', {'auto', 'spin', 'msr', 'permutation'}, ...
                                'Value', 'auto');
        S.nullDrop.Layout.Row = 4;
        S.nullDrop.Layout.Column = 2;
        S.h.labPermV = uilabel(gv, 'Text', T('perms'), 'HorizontalAlignment', 'right');
        S.h.labPermV.Layout.Row = 5;
        S.permSpinner = uispinner(gv, 'Limits', [20 5000], 'Step', 50, 'Value', 500);
        S.permSpinner.Layout.Row = 5;
        S.permSpinner.Layout.Column = 2;
        S.h.varRunBtn = uibutton(gv, 'push', 'FontWeight', 'bold', ...
            'BackgroundColor', [0.16 0.42 0.36], 'FontColor', [1 1 1], ...
            'ButtonPushedFcn', @onVarRun);
        S.h.varRunBtn.Layout.Row = 7;
        S.h.varRunBtn.Layout.Column = [1 2];

        resTabs = uitabgroup(g);
        resTabs.Layout.Column = 2;
        S.h.tAxis = uitab(resTabs);
        ga = uigridlayout(S.h.tAxis, [1 1]);
        S.axisTable = uitable(ga, 'ColumnName', 'auto', 'RowName', []);

        S.h.tProf = uitab(resTabs);
        gp1 = uigridlayout(S.h.tProf, [1 1]);
        S.varAxes = uiaxes(gp1);
        styleAxes(S.varAxes);

        S.h.tLayer = uitab(resTabs);
        gl = uigridlayout(S.h.tLayer, [1 2]);
        S.layerAxes = uiaxes(gl);
        styleAxes(S.layerAxes);
        S.epiAxes = uiaxes(gl);
        styleAxes(S.epiAxes);

        S.h.tSweep = uitab(resTabs);
        gs2 = uigridlayout(S.h.tSweep, [1 1]);
        S.sweepAxes = uiaxes(gs2);
        styleAxes(S.sweepAxes);

        S.h.resTabs = resTabs;
    end

    function onVarRun(~, ~)
        S.h.varRunBtn.Enable = 'off';
        S.progress.Text = T('preparing');
        drawnow;
        try
            requireCohort();
            if S.coh.n < 2
                uialert(fig, sprintf(T('oneSubjectBlock'), S.coh.n), T('analyzeFail'));
                return
            end
            metrics = S.varMetrics.Value;
            if ischar(metrics) || isstring(metrics)
                metrics = cellstr(metrics);
            end
            if isempty(metrics)
                uialert(fig, T('pickOne'), T('analyzeFail'));
                return
            end
            est = S.estDrop.Value;
            if strcmp(est, 'all')
                estList = {'regional', 'node', 'edge'};
            else
                estList = {est};
            end
            bundle = sfcx.analyze(S.coh, metrics, 'estimators', estList, ...
                                  'nPerm', S.permSpinner.Value, ...
                                  'null', S.nullDrop.Value, 'seed', 0, ...
                                  'progress', @onProgress);
            S.bundle = bundle;
            refreshVarViews();
            S.h.resTabs.SelectedTab = S.h.tAxis;
            S.progress.Text = sprintf(T('doneFmt'), numel(metrics), ...
                numel(estList), bundle.runtimeS);
            logLine(sprintf('variability: %d metrics, %.1f s', ...
                            numel(metrics), bundle.runtimeS));
        catch err
            S.progress.Text = T('failed');
            uialert(fig, err.message, T('analyzeFail'));
        end
        S.h.varRunBtn.Enable = 'on';
    end

    function refreshVarViews()
        b = S.bundle;
        rows = cell(0, 6);
        for i = 1:numel(b.metrics)
            for j = 1:numel(b.metrics(i).estimators)
                e = b.metrics(i).estimators(j);
                rows(end + 1, :) = {b.metrics(i).name, e.name, ...
                                    e.axis.r, e.axis.p, e.axis.n, e.axis.null}; %#ok<AGROW>
            end
        end
        S.axisTable.Data = rows;

        cla(S.varAxes); hold(S.varAxes, 'on');
        legendNames = cell(1, numel(b.metrics));
        for i = 1:numel(b.metrics)
            idx = strcmp({b.metrics(i).estimators.name}, b.metrics(i).primaryEstimator);
            v = b.metrics(i).estimators(idx).values;
            plot(S.varAxes, 1:numel(v), v, 'LineWidth', 1.2);
            legendNames{i} = b.metrics(i).name;
        end
        hold(S.varAxes, 'off');
        grid(S.varAxes, 'on');
        legend(S.varAxes, legendNames, 'Location', 'best');
        xlabel(S.varAxes, T('profileX'));
        ylabel(S.varAxes, T('profileY'));
        title(S.varAxes, T('profileTitle'));

        cla(S.layerAxes);
        vals = zeros(numel(b.metrics), 2);
        names = cell(1, numel(b.metrics));
        for i = 1:numel(b.metrics)
            vals(i, :) = [b.metrics(i).layers.within, b.metrics(i).layers.between];
            names{i} = b.metrics(i).name;
        end
        bar(S.layerAxes, vals);
        S.layerAxes.XTick = 1:numel(names);
        S.layerAxes.XTickLabel = names;
        S.layerAxes.XTickLabelRotation = 30;
        legend(S.layerAxes, {T('within'), T('between')}, 'Location', 'best');
        grid(S.layerAxes, 'on');
        title(S.layerAxes, T('layerTitle'));

        cla(S.epiAxes); hold(S.epiAxes, 'on');
        i0 = 1;
        e0 = b.metrics(i0).estimators(strcmp({b.metrics(i0).estimators.name}, ...
                                             b.metrics(i0).primaryEstimator));
        v = e0.values;
        plot(S.epiAxes, 1:numel(v), v, 'Color', [0.45 0.55 0.62], 'LineWidth', 1);
        focal = b.metrics(i0).epicentres.regions;
        if ~isempty(focal)
            scatter(S.epiAxes, focal, v(focal), 46, 'filled', ...
                    'MarkerFaceColor', [0.76 0.30 0.36]);
        end
        hold(S.epiAxes, 'off');
        grid(S.epiAxes, 'on');
        title(S.epiAxes, sprintf(T('epiTitleFmt'), numel(focal), b.metrics(i0).name));

        cla(S.sweepAxes);
        imagesc(S.sweepAxes, b.sweep.r);
        colormap(S.sweepAxes, sfcxColors());
        colorbar(S.sweepAxes);
        S.sweepAxes.XTick = 1:numel(b.sweep.cols);
        S.sweepAxes.XTickLabel = b.sweep.cols;
        S.sweepAxes.YTick = 1:numel(b.sweep.rows);
        S.sweepAxes.YTickLabel = b.sweep.rows;
        clim(S.sweepAxes, [-1 1]);
        title(S.sweepAxes, T('sweepTitle'));
    end

% ===================================================================
%  TAB 3  COUPLING VALUES
% ===================================================================
    function buildMagTab()
        g = uigridlayout(S.h.tMag, [1 2]);
        g.ColumnWidth = {320, '1x'};
        g.ColumnSpacing = 12;
        g.Padding = [6 6 6 6];

        S.h.pMagOpts = uipanel(g, 'Title', T('secMagOpts'), 'BackgroundColor', [1 1 1]);
        S.h.pMagOpts.Layout.Column = 1;
        gm2 = uigridlayout(S.h.pMagOpts, [6 2]);
        gm2.ColumnWidth = {86, '1x'};
        gm2.RowHeight = {18, 200, 24, 24, 12, 40};
        gm2.RowSpacing = 6;
        gm2.Padding = [8 8 8 8];

        S.h.magMetricHint = uilabel(gm2, 'Text', T('magMetricHint'), ...
            'FontColor', [0.42 0.49 0.55]);
        S.h.magMetricHint.Layout.Row = 1;
        S.h.magMetricHint.Layout.Column = [1 2];
        S.magMetrics = uilistbox(gm2, 'Multiselect', 'on');
        S.magMetrics.Layout.Row = 2;
        S.magMetrics.Layout.Column = [1 2];
        S.h.magCsv = uicheckbox(gm2, 'Text', T('csvLab'), 'Value', true);
        S.h.magCsv.Layout.Row = 3;
        S.h.magCsv.Layout.Column = [1 2];
        S.h.magMat = uicheckbox(gm2, 'Text', T('matLab'), 'Value', true);
        S.h.magMat.Layout.Row = 4;
        S.h.magMat.Layout.Column = [1 2];
        runRow = uigridlayout(gm2);
        runRow.Layout.Row = 6;
        runRow.Layout.Column = [1 2];
        runRow.ColumnWidth = {'1x', '1x'};
        runRow.RowHeight = {40};
        runRow.Padding = [0 0 0 0];
        S.h.magRunBtn = uibutton(runRow, 'push', 'FontWeight', 'bold', ...
            'BackgroundColor', [0.16 0.42 0.36], 'FontColor', [1 1 1], ...
            'ButtonPushedFcn', @onMagRun);
        S.h.magRunBtn.Layout.Column = 1;
        S.h.magExportBtn = uibutton(runRow, 'push', 'ButtonPushedFcn', @onMagExport);
        S.h.magExportBtn.Layout.Column = 2;

        magTabs = uitabgroup(g);
        magTabs.Layout.Column = 2;
        S.h.tMeanMat = uitab(magTabs);
        gm3 = uigridlayout(S.h.tMeanMat, [2 1]);
        gm3.RowHeight = {'1x', 26};
        S.matAxes = uiaxes(gm3);
        styleAxes(S.matAxes);
        S.matInfo = uilabel(gm3, 'FontColor', [0.42 0.49 0.55]);

        S.h.tRegional = uitab(magTabs);
        gr = uigridlayout(S.h.tRegional, [1 1]);
        S.regAxes = uiaxes(gr);
        styleAxes(S.regAxes);

        S.h.magTabs = magTabs;
    end

    function onMagRun(~, ~)
        S.h.magRunBtn.Enable = 'off';
        S.progress.Text = T('preparing');
        drawnow;
        try
            requireCohort();
            metrics = S.magMetrics.Value;
            if ischar(metrics) || isstring(metrics)
                metrics = cellstr(metrics);
            end
            if isempty(metrics)
                uialert(fig, T('pickOne'), T('magFail'));
                return
            end
            vals = sfcx.couplingValues(S.coh, metrics);
            S.vals = vals;
            refreshMagViews();
            S.h.magTabs.SelectedTab = S.h.tMeanMat;
            S.progress.Text = sprintf(T('magDoneFmt'), numel(metrics));
            logLine(sprintf('coupling values: %s (%d subjects)', ...
                strjoin(metrics, ', '), S.coh.n));
        catch err
            S.progress.Text = T('failed');
            uialert(fig, err.message, T('magFail'));
        end
        S.h.magRunBtn.Enable = 'on';
    end

    function refreshMagViews()
        cla(S.matAxes);
        matNames = {};
        for i = 1:numel(S.vals)
            if ~isempty(S.vals(i).meanMatrix)
                matNames{end + 1} = S.vals(i).name; %#ok<AGROW>
            end
        end
        if isempty(matNames)
            S.matInfo.Text = T('noMatrix');
            imagesc(S.matAxes, nan(10));
        else
            i0 = find(strcmp({S.vals.name}, matNames{1}), 1);
            M = S.vals(i0).meanMatrix;
            imagesc(S.matAxes, M);
            colormap(S.matAxes, sfcxColors());
            colorbar(S.matAxes);
            axis(S.matAxes, 'square');
            title(S.matAxes, sprintf(T('meanMatTitle'), S.vals(i0).name));
            S.matInfo.Text = sprintf('%s: %s', T('tabMeanMat'), strjoin(matNames, ', '));
        end

        cla(S.regAxes);
        R = [];
        counts = [];
        for i = 1:numel(S.vals)
            R = [R; S.vals(i).regional]; %#ok<AGROW>
            counts(end + 1) = size(S.vals(i).regional, 1); %#ok<AGROW>
        end
        imagesc(S.regAxes, R);
        colormap(S.regAxes, parula);
        colorbar(S.regAxes);
        S.regAxes.YTick = 1:size(R, 1);
        labels = strcat(repelem({S.vals.name}, counts));
        S.regAxes.YTickLabel = labels;
        if S.coh.n == 1
            xlabel(S.regAxes, T('regionalXSingle'));
        else
            xlabel(S.regAxes, T('regionalX'));
        end
        ylabel(S.regAxes, T('regionalY'));
        title(S.regAxes, T('regionalTitle'));
    end

    function onMagExport(~, ~)
        if isempty(S.vals)
            uialert(fig, T('magRunFirst'), T('export'));
            return
        end
        d = uigetdir(pwd, T('pickDir'));
        if isequal(d, 0)
            return
        end
        fmt = {};
        if S.h.magCsv.Value
            fmt{end + 1} = 'csv';
        end
        if S.h.magMat.Value
            fmt{end + 1} = 'mat';
        end
        if isempty(fmt)
            uialert(fig, T('pickFormat'), T('export'));
            return
        end
        try
            sfcx.exportCouplingValues(S.vals, d, fmt);
            S.progress.Text = sprintf(T('exportedFmt'), d);
            logLine(sprintf('exported coupling values -> %s', d));
        catch err
            uialert(fig, err.message, T('export'));
        end
    end

% ===================================================================
%  TAB 4  LOG
% ===================================================================
    function buildLogTab()
        g = uigridlayout(S.h.tLog, [2 1]);
        g.RowHeight = {'1x', 44};
        g.Padding = [6 6 6 6];
        S.logArea = uitextarea(g, 'Editable', 'off', 'FontName', 'Consolas');
        S.h.exportVarBtn = uibutton(g, 'push', 'Text', T('exportVar'), ...
            'ButtonPushedFcn', @onExportVar);
    end

    function onExportVar(~, ~)
        if isempty(S.bundle)
            uialert(fig, T('runFirst'), T('export'));
            return
        end
        d = uigetdir(pwd, T('pickDir'));
        if isequal(d, 0)
            return
        end
        sfcx.exportResults(S.bundle, d, S.coh);
        S.progress.Text = sprintf(T('exportedFmt'), d);
        logLine(sprintf('exported variability results -> %s', d));
    end

    function logLine(txt)
        txt = char(txt);
        old = S.logArea.Value;
        if ischar(old)
            old = cellstr(old);
        end
        old = old(:);
        S.logArea.Value = [old; {txt}]; %#ok<AGROW>
    end

% ------------------------------- shared ----------------------------
    function requireCohort()
        if isempty(S.coh)
            error('sfcx:app', '%s', T('noCohort'));
        end
    end

    function onProgress(frac, msg)
        S.progress.Text = sprintf('%3.0f%%  %s', frac * 100, msg);
        drawnow limitrate;
    end

    function styleAxes(ax)
        ax.Box = 'off';
        ax.FontName = 'Helvetica';
        ax.FontSize = 9;
    end

    function applyLang()
        fig.Name = T('appName');
        S.h.rb1.Text = T('modeSingle');
        S.h.rb2.Text = T('modeMulti');
        S.h.langBtn.Text = T('lang');
        S.h.tData.Title = T('tabData');
        S.h.tVar.Title = T('tabVar');
        S.h.tMag.Title = T('tabMag');
        S.h.tLog.Title = T('tabLog');
        S.h.pSlot.Title = T('secSingle');
        S.h.pBuild.Title = T('secBuild');
        S.h.pVarOpts.Title = T('secVarOpts');
        S.h.pMagOpts.Title = T('secMagOpts');
        S.h.demoBtn.Text = T('demo');
        S.h.annoBtn.Text = T('annoLab');
        S.h.buildBtn.Text = T('build');
        S.h.varRunBtn.Text = T('runVar');
        S.h.magRunBtn.Text = T('runMag');
        S.h.magExportBtn.Text = T('export');
        S.h.exportVarBtn.Text = T('exportVar');
        S.h.labEstV.Text = T('estimator');
        S.h.labNullV.Text = T('nullModel');
        S.h.labPermV.Text = T('perms');
        S.h.varMetricHint.Text = T('metricHint');
        S.h.magMetricHint.Text = T('magMetricHint');
        S.h.magCsv.Text = T('csvLab');
        S.h.magMat.Text = T('matLab');
        S.h.tAxis.Title = T('tabAxis');
        S.h.tProf.Title = T('tabProfile');
        S.h.tLayer.Title = T('tabLayer');
        S.h.tSweep.Title = T('tabSweep');
        S.h.tMeanMat.Title = T('tabMeanMat');
        S.h.tRegional.Title = T('tabRegional');
        S.axisTable.ColumnName = T('axisCols');
        S.h.slotLabsc.Text = T('scLab');
        S.h.slotLabfc.Text = T('fcLab');
        S.h.slotLabbold.Text = T('boldLab');
        S.h.slotLabscBatch.Text = T('scBatchLab');
        S.h.slotLabfcBatch.Text = T('fcBatchLab');
        S.h.slotLabboldBatch.Text = T('boldBatchLab');
        S.h.slotLabnetworks.Text = T('netLab');
        S.h.slotLabhierarchy.Text = T('hierLab');
        S.h.slotLabcoords.Text = T('coordLab');
        S.h.slotLabbehaviour.Text = T('behLab');
        S.h.slotLabcovariates.Text = T('covLab');
        updateModePanels();
        refreshMetricList(S.varMetrics);
        refreshMetricList(S.magMetrics);
    end

    function refreshMetricList(box)
        cat = sfcx.catalog();
        ids = {cat.metrics.id};
        disp_ = cell(size(ids));
        for i = 1:numel(ids)
            if strcmp(S.lang, 'zh') && isfield(cat.metrics, 'label_zh') && ...
                    ~isempty(cat.metrics(i).label_zh)
                disp_{i} = sprintf('%s  (%s)', cat.metrics(i).label_zh, ids{i});
            else
                disp_{i} = sprintf('%s  (%s)', cat.metrics(i).label, ids{i});
            end
        end
        keep = box.Value;
        box.Items = disp_;
        box.ItemsData = ids;
        if isempty(keep) || isempty(intersect(keep, ids, 'stable'))
            box.Value = ids(1);
        else
            box.Value = intersect(keep, ids, 'stable');
        end
    end

    function v = readOne(full, kind)
        [~, ~, ext] = fileparts(full);
        if strcmpi(ext, '.mat')
            s = load(full);
            v = fieldPick(s, kind);
        else
            v = readmatrix(full);
        end
        v = double(v);
        if strcmp(kind, 'bold') && ismatrix(v) && size(v, 1) > size(v, 2)
            v = v';  % time x nodes -> nodes x time
        end
    end

    function v = loadAny(prompt)
        v = [];
        [file, path] = uigetfile({'*.mat', 'MAT'}, prompt);
        if isequal(file, 0)
            return
        end
        v = load(fullfile(path, file));
    end

    function v = matTo3D(full)
        s = load(full);
        v = [];
        fn = fieldnames(s);
        for i = 1:numel(fn)
            x = s.(fn{i});
            if iscell(x)
                first = double(x{1});
                v = zeros([numel(x), size(first)]);
                for k = 1:numel(x)
                    v(k, :, :, :) = double(x{k});
                end
                return
            elseif isnumeric(x) && ndims(x) == 3
                v = double(x);
                return
            end
        end
        if isempty(v)
            error('sfcx:app', 'no 3-D numeric array or cell array in %s', full);
        end
    end

    function v = fieldPick(s, kind)
        v = [];
        prefs = {'sc', 'SC', 'struct', 'structural'; ...
                 'fc', 'FC', 'func', 'functional'; ...
                 'bold', 'BOLD', 'timeseries', 'ts'};
        row = strcmp(prefs(:, 1), kind);
        aliases = prefs{find(row, 1), :};
        for i = 1:numel(aliases)
            if isfield(s, aliases{i})
                v = s.(aliases{i});
                return
            end
        end
        fn = fieldnames(s);
        if isscalar(fn)
            v = s.(fn{1});
        else
            error('sfcx:app', 'no %s field; available: %s', kind, strjoin(fn', ', '));
        end
    end

    function txt = sizeText(v)
        dims = size(v);
        txt = strjoin(arrayfun(@(d) sprintf('%d', d), dims, ...
                               'UniformOutput', false), '×');
    end

    function s = T(key)
        s = strings_(S.lang, key);
    end
end

% ======================================================================
%  helpers outside the UI closure
% ======================================================================
function cmap = sfcxColors()
cmap = [linspace(0.13, 0.99, 256)', linspace(0.35, 0.99, 256)', linspace(0.42, 0.99, 256)'];
end

function s = strings_(lang, key)
%STRINGS_ UI strings, zh and en.
persistent ZH EN
if isempty(ZH)
    ZH = containers.Map();
    EN = containers.Map();
    pair('appName', 'SFCX  |  结构–功能耦合工具箱', 'SFCX  |  Structure–function coupling toolbox');
    pair('tabData', '① 数据导入', '① Data');
    pair('tabVar', '② 变异分析', '② Variability');
    pair('tabMag', '③ 耦合值计算', '③ Coupling values');
    pair('tabLog', '④ 日志', '④ Log');
    pair('modeSingle', '单个被试模式', 'Single subject');
    pair('modeMulti', '多被试批量模式', 'Multi-subject batch');
    pair('modeSingleHint', '导入一名被试的 SC + FC/BOLD。耦合值计算（页签③）可用；变异分析（页签②）需要 ≥2 名被试。', ...
         'One subject''s SC + FC/BOLD. Coupling values (tab ③) work; variability (tab ②) needs ≥2 subjects.');
    pair('modeMultiHint', '选择一个文件夹（每人一个文件）或一个含三维数组的 .mat；支持 SC、FC 或 BOLD。', ...
         'Pick a folder (one file per subject) or a .mat with 3-D arrays; SC, FC or BOLD.');
    pair('secSingle', '单个被试数据', 'Single-subject data');
    pair('secBatch', '多被试批量数据', 'Multi-subject batch data');
    pair('secBuild', '注释（可选）', 'Annotation (optional)');
    pair('scLab', 'SC 矩阵', 'SC matrix');
    pair('fcLab', 'FC 矩阵', 'FC matrix');
    pair('boldLab', 'BOLD 序列', 'BOLD series');
    pair('scBatchLab', 'SC 批量', 'SC batch');
    pair('fcBatchLab', 'FC 批量', 'FC batch');
    pair('boldBatchLab', 'BOLD 批量', 'BOLD batch');
    pair('netLab', '网络标签', 'Networks');
    pair('hierLab', '层级', 'Hierarchy');
    pair('coordLab', '坐标', 'Coords');
    pair('behLab', '行为', 'Behaviour');
    pair('covLab', '协变量', 'Covariates');
    pair('annoLab', '注释 .mat（全部字段）', 'Annotation .mat (all fields)');
    pair('load', '载入…', 'Load…');
    pair('empty', '—', '—');
    pair('demo', '示例数据', 'Demo data');
    pair('build', '构建队列', 'Build cohort');
    pair('builtFmt', '队列就绪：%d 名被试 × %d 脑区（%d 条边）', ...
         'Cohort ready: %d subjects × %d nodes (%d edges)');
    pair('oneSubjectWarn', '只有 %d 名被试：可计算耦合值（页签③），变异分析（页签②）不可用。', ...
         'Only %d subject(s): coupling values (tab ③) work, variability (tab ②) does not.');
    pair('smallWarn', '被试数 %d < 8：变异统计噪声较大，结果仅作演示参考。', ...
         'n = %d < 8: variability estimates are noisy; treat results as a demo.');
    pair('oneSubjectBlock', '当前队列只有 %d 名被试，变异度统计在数学上未定义。请使用页签③计算耦合值，或导入多被试数据。', ...
         'This cohort has %d subject(s); dispersion is undefined. Use tab ③ for coupling values or import more subjects.');
    pair('needData', '请先导入 SC / FC / BOLD（或点击示例数据）', ...
         'Import SC / FC / BOLD first (or load the demo)');
    pair('buildFail', '构建失败', 'Build failed');
    pair('batchChoice', '批量导入方式？', 'Batch import method?');
    pair('batchTitle', '批量导入', 'Batch import');
    pair('batchFolder', '文件夹（每人一个文件）', 'Folder (one file per subject)');
    pair('batchMat', '单个 .mat（三维数组/元胞）', 'One .mat (3-D array / cell)');
    pair('cancel', '取消', 'Cancel');
    pair('pickBatchDir', '选择被试文件夹', 'Choose the subjects folder');
    pair('pickBatchMat', '选择 .mat 文件', 'Choose a .mat file');
    pair('pickSC', '选择结构连接矩阵', 'Choose the SC matrix');
    pair('pickFC', '选择功能连接矩阵', 'Choose the FC matrix');
    pair('pickBold', '选择 BOLD 时间序列', 'Choose the BOLD series');
    pair('pickAnno', '选择注释文件', 'Choose the annotation file');
    pair('annoNoneLong', '文件中没有 networks / hierarchy / coords / behaviour / covariates', ...
         'The file has none of networks / hierarchy / coords / behaviour / covariates');
    pair('smartFilled', '已导入注释：%s', 'Imported annotation: %s');
    pair('loadTitle', '载入', 'Load');
    pair('demoReady', '示例数据就绪（160 名被试 × 100 脑区）', ...
         'Demo data ready (160 subjects × 100 regions)');
    pair('secVarOpts', '变异分析设置（需 ≥2 名被试）', 'Variability options (n ≥ 2)');
    pair('metricHint', '耦合定义（可多选）', 'Coupling definitions (multi-select)');
    pair('estimator', '估计量', 'Estimator');
    pair('nullModel', '零模型', 'Null model');
    pair('perms', '置换次数', 'Permutations');
    pair('runVar', '运行变异分析', 'Run variability');
    pair('tabAxis', '连接轴', 'Axis');
    pair('tabProfile', '变异剖面', 'Profile');
    pair('tabLayer', '层级与震中', 'Layers & epicentres');
    pair('tabSweep', '坐标扫描', 'Sweep');
    pair('profileTitle', '各定义的离散度剖面', 'Dispersion profile');
    pair('profileX', '区域索引', 'Region index');
    pair('profileY', '离散度（SD）', 'Dispersion (SD)');
    pair('layerTitle', '网络内 vs 网络间', 'Within vs between');
    pair('within', '网络内', 'within');
    pair('between', '网络间', 'between');
    pair('epiTitleFmt', '震中 %d 个（%s）', '%d epicentres (%s)');
    pair('sweepTitle', '定义 × 估计量的轴相关系数', 'Axis r: definition × estimator');
    pair('axisCols', {'定义', '估计量', 'r', 'p', 'n', '零模型'}, ...
         {'metric', 'estimator', 'r', 'p', 'n', 'null'});
    pair('secMagOpts', '耦合值计算（单/多被试均可）', 'Coupling values (single or batch)');
    pair('magMetricHint', '耦合定义（可多选）', 'Coupling definitions (multi-select)');
    pair('csvLab', '导出 CSV（逐被试值）', 'Export CSV (per-subject values)');
    pair('matLab', '导出 MAT（含平均矩阵）', 'Export MAT (with mean matrix)');
    pair('runMag', '计算耦合值', 'Compute values');
    pair('tabMeanMat', '平均耦合矩阵', 'Mean matrix');
    pair('tabRegional', '区域耦合热图', 'Regional heat map');
    pair('meanMatTitle', '组平均耦合矩阵（%s）', 'Mean coupling matrix (%s)');
    pair('noMatrix', '所选定义均为区域级输出（如 Granger），没有 P×P 矩阵；请看「区域耦合热图」。', ...
         'All selected definitions are regional (e.g. Granger): no P×P matrix; see the regional heat map.');
    pair('regionalTitle', '各定义的区域级耦合值', 'Regional coupling per definition');
    pair('regionalX', '脑区', 'Region');
    pair('regionalXSingle', '脑区（单被试）', 'Region (single subject)');
    pair('regionalY', '定义 / 被试', 'Definition / subject');
    pair('magDoneFmt', '耦合值计算完成：%d 个定义。可点击「导出结果」保存。', ...
         'Coupling values computed: %d definitions. Press Export to save.');
    pair('magRunFirst', '先点击「计算耦合值」', 'Press Compute values first');
    pair('magFail', '耦合值计算失败', 'Coupling computation failed');
    pair('pickFormat', '请至少选择一种导出格式', 'Select at least one export format');
    pair('export', '导出结果', 'Export results');
    pair('exportVar', '导出变异分析结果（CSV/MAT）', 'Export variability results (CSV/MAT)');
    pair('pickDir', '选择导出目录', 'Choose export folder');
    pair('exportedFmt', '已导出到 %s', 'Exported to %s');
    pair('runFirst', '先运行一次分析', 'Run an analysis first');
    pair('noCohort', '请先在「数据导入」页签构建队列', 'Build a cohort in the Data tab first');
    pair('pickOne', '请至少选择一个耦合定义', 'Select at least one definition');
    pair('preparing', '准备数据…', 'Preparing…');
    pair('doneFmt', '完成：%d 个定义 × %d 个估计量，%.1f 秒', ...
         'Done: %d definitions × %d estimators, %.1f s');
    pair('failed', '失败', 'Failed');
    pair('analyzeFail', '分析失败', 'Analysis failed');
    pair('ready', '就绪 — 先导入数据并构建队列', 'Ready — import data and build a cohort');
    pair('lang', 'English', '中文');
    pair('logHead', '提示：单被试只能用页签③；变异分析需 ≥2 名被试。', ...
         'Note: single subjects work in tab ③ only; variability needs ≥2 subjects.');
end
if strcmp(lang, 'en') && EN.isKey(key)
    s = EN(key);
elseif ZH.isKey(key)
    s = ZH(key);
else
    s = key;
end

    function pair(k, zh, en)
        ZH(k) = zh;
        EN(k) = en;
    end
end
