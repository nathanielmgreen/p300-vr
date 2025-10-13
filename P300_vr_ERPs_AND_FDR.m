clc; clear; close all;

%% =============================
%            SETUP
% ==============================
mat_file    = 'epochs_data_08302025_new.mat';   % <-- path to your .mat
chan_labels = {'Cz','Pz','Oz','O1','O2'};  % columns in each row of the 3x5 grid
baseline_ms = [-300 0];                    % baseline window (ms); set [] to skip
xlim_ms     = [-300 1000];                 % plot window
fdr_win     = [200 600];                     % FDR correction window (e.g., P300-focused)
alpha_raw   = 0.05;                        % raw p threshold
q_fdr       = 0.05;                        % FDR (BH) threshold
stats_mode  = 'trial_level';               % 'trial_level' or 'paired_by_file' (no file IDs here, so use trial_level)

% --- Trimming policy knobs to match your onset decisions ---
INNER_FIX_MAX_MS  = 400;  % inner: cut fixation after 400 ms (tune to 450, etc.)
OUTER_TAIL_MAX_MS = 500;  % outer: cut BOTH GIP and FIX after 500 ms
% GIP/Fix min windows (keep your previous values)
rules_inner.gip = struct('min_ms', 50,  'max_ms', 400);
rules_inner.fix = struct('min_ms', 80,  'max_ms', 700);
rules_outer.gip = struct('min_ms', 50,  'max_ms', 600);
rules_outer.fix = struct('min_ms', 120, 'max_ms', 1200);
% Apply your new tail narrowing:
rules_inner.fix.max_ms = INNER_FIX_MAX_MS;
rules_outer.gip.max_ms = min(rules_outer.gip.max_ms, OUTER_TAIL_MAX_MS);
rules_outer.fix.max_ms = min(rules_outer.fix.max_ms, OUTER_TAIL_MAX_MS);

% Which locks use rules when selecting trials:
apply_rules = struct('gip', false, 'fix', false, 'stim', false); % Stimulus stays rule-free

%% =============================
%        LOAD EPOCHS STRUCT
% ==============================
S = load(mat_file); 
if ~isfield(S, 'epochs'), error('File does not contain variable ''epochs''.'); end
epochs = S.epochs;

% Expect fields: epochs.stimulus_onset, epochs.gip_onset, epochs.fix_onset
% Each has inner_data (ch x time x trials), inner_stimuli (1 x trials cell), inner_timeAfterStim (1 x trials),
% and outer_* counterparts. Also epochs.times (1 x T), epochs.chanlocs (1 x N struct with .labels)

times = epochs.times(:)';  % 1 x T (likely -500..1800 ms)

% Map requested channel labels -> indices in chanlocs
chIdxMap = map_channel_indices(epochs.chanlocs, chan_labels);

% --- Sanity check: inner vs outer for GIP onset ---
X_in  = epochs.gip_onset.inner_data;   % ch x time x trials
X_out = epochs.gip_onset.outer_data;   % ch x time x trials

fprintf('GIP inner size:  %s\n', mat2str(size(X_in)));
fprintf('GIP outer size:  %s\n', mat2str(size(X_out)));

% Compare a stable summary (mean over trials, then over time) per channel
m_in  = squeeze(mean(mean(X_in,  2, 'omitnan'), 3, 'omitnan')); % ch x 1
m_out = squeeze(mean(mean(X_out, 2, 'omitnan'), 3, 'omitnan')); % ch x 1

same_size = isequal(size(X_in), size(X_out));
same_vals = same_size && all(abs(m_in - m_out) < 1e-10 | (isnan(m_in) & isnan(m_out)));

fprintf('Inner/Outer GIP means identical? %d (1=true)\n', same_vals);

% Optional: correlation of channel-wise means
if same_size
    R = corr(m_in, m_out, 'rows','pairwise');
    fprintf('Corr(inner,outer) over channels: %.4f\n', R);
end


%% =============================
%      BUILD ERPs (INNER/OUTER)
% ==============================
ERP_INNER = build_erps_from_epochs(epochs, 'inner', times, chIdxMap, chan_labels, baseline_ms, rules_inner, apply_rules);
ERP_OUTER = build_erps_from_epochs(epochs, 'outer', times, chIdxMap, chan_labels, baseline_ms, rules_outer, apply_rules);

%% =============================
%        PLOTTING (BOTH)
% ==============================
% Plain figures (no stats)
plot_erps_std_dev_multichan(ERP_INNER, 'Eye Shifting (Inner)', chan_labels, xlim_ms);
plot_erps_std_dev_multichan(ERP_OUTER, 'Head Turning (Outer)', chan_labels, xlim_ms);

% With stats overlays (gray=raw p<.05, black=FDR q<.05) within fdr_win
plot_erps_std_dev_multichan_with_stats(ERP_INNER, 'Eye Shifting (Inner)', chan_labels, xlim_ms, alpha_raw, q_fdr, stats_mode, fdr_win);
plot_erps_std_dev_multichan_with_stats(ERP_OUTER, 'Head Turning (Outer)', chan_labels, xlim_ms, alpha_raw, q_fdr, stats_mode, fdr_win);

disp('Done.');

%% =============================
%            FUNCTIONS
% ==============================

function ERP = build_erps_from_epochs(epochs, side, times, chIdxMap, chan_labels, baseline_ms, rules, apply_rules)
% side: 'inner' or 'outer'
% Builds trial-level and summary ERPs for locks: stimulus_onset (no rules), gip_onset (rules), fix_onset (rules)

locks_src = {'gip_onset','fix_onset','stimulus_onset'};
locks_out = {'gip','fix','stim'};
ERP = struct();

for L = 1:numel(locks_src)
    src = locks_src{L};
    out = locks_out{L};

    % Pull data arrays from struct: data: ch x time x trials
    data   = epochs.(src).([side '_data']);
    labels = epochs.(src).([side '_stimuli']);       % 1 x trials cell ('standard'/'deviant')
    tAfter = epochs.(src).([side '_timeAfterStim']); % 1 x trials (ms) (0 for stimulus_onset)

    % DEBUG: show counts before/after trimming for this lock
    use_rules = apply_rules.(out);
    if isfield(rules, out), this_rule = rules.(out); else, this_rule = struct(); end
    debug_print_counts(src, side, labels, tAfter, this_rule, use_rules);

    % --- Trial selection by condition ---
    is_std = strcmpi(labels, 'standard');
    is_dev = strcmpi(labels, 'deviant');

    % --- Apply trimming rules (only for gip/fix if apply_rules says so) ---
    use_rules = apply_rules.(out); % out is 'gip'|'fix'|'stim'
    if use_rules
        rule = rules.(out);
        keep = true(size(tAfter));
        if isfield(rule,'min_ms') && ~isempty(rule.min_ms), keep = keep & (tAfter >= rule.min_ms); end
        if isfield(rule,'max_ms') && ~isempty(rule.max_ms), keep = keep & (tAfter <= rule.max_ms); end
        is_std = is_std & keep;
        is_dev = is_dev & keep;
    end
    % For stimulus-locked, use_rules is false per your policy → no trimming

    % Optional: baseline-correct (per channel, per trial)
    if ~isempty(baseline_ms)
        data = baseline_epochs(data, times, baseline_ms);
    end

    % Initialize stores
    ERP.(out).times = times;
    for c = 1:numel(chan_labels)
        lab = chan_labels{c};
        idx = chIdxMap.(lab);
        if isempty(idx)
            warning('Channel %s not found; leaving empty.', lab);
            ERP.trials.(out).std.(lab) = [];
            ERP.trials.(out).dev.(lab) = [];
            ERP.(out).std.(lab) = summarize_trials_struct([]);
            ERP.(out).dev.(lab) = summarize_trials_struct([]);
            continue;
        end

        % Extract trials x time
        X = squeeze(data(idx, :, :))';  % trials x time

        X_std = X(is_std, :);
        X_dev = X(is_dev, :);

        % Save trials for stats
        ERP.trials.(out).std.(lab) = X_std;
        ERP.trials.(out).dev.(lab) = X_dev;

        % Summaries
        ERP.(out).std.(lab) = summarize_trials_struct(X_std);
        ERP.(out).dev.(lab) = summarize_trials_struct(X_dev);
    end
end
end

function chIdxMap = map_channel_indices(chanlocs, wanted_labels)
% Map wanted EEG labels to indices using chanlocs struct (1xN with .labels)
chIdxMap = struct();
have_labels = arrayfun(@(x) string(x.labels), chanlocs, 'UniformOutput', true);
for k = 1:numel(wanted_labels)
    lab = wanted_labels{k};
    idx = find(strcmpi(have_labels, lab), 1);
    if isempty(idx)
        warning('Channel %s not found in chanlocs.', lab);
        chIdxMap.(lab) = [];
    else
        chIdxMap.(lab) = idx;
    end
end
end

function data_bc = baseline_epochs(data, times, baseline_ms)
% Baseline-correct epochs: data (ch x time x trials), subtract mean over baseline_ms per channel/trial.
data_bc = data;
t = times(:)';
bmin = max(baseline_ms(1), t(1));
bmax = min(baseline_ms(2), t(end));
if ~(bmin < bmax)
    warning('baseline_epochs: invalid baseline window; skipping baseline.');
    return;
end
bmask = (t >= bmin & t <= bmax);
if ~any(bmask)
    warning('baseline_epochs: no baseline samples; skipping baseline.');
    return;
end
for tr = 1:size(data,3)
    bl = mean(data(:, bmask, tr), 2, 'omitnan');  % ch x 1
    data_bc(:, :, tr) = data(:, :, tr) - bl;
end
end

function S = summarize_trials_struct(X)
% trials x time -> mean/SEM/n (handles empty)
if isempty(X)
    S.mean = []; S.sem = []; S.n = 0;
else
    S.mean = mean(X,1, 'omitnan');
    S.sem  = std(X,0,1, 'omitnan') / sqrt(size(X,1));
    S.n    = size(X,1);
end
end

% ======== PLOTTERS (reuse your existing ones) ========

function plot_erps_std_dev_multichan(ERP, panel_title, chan_labels, xlim_ms)
% 3 rows (GIP/Fix/Stim), columns = chan_labels; clean legend; bottom-row xlabels only.
ylims = [-4 4]; % e.g., [-4 4]; leave [] to autoscale
locks = {'gip','fix','stim'}; row_titles = {'GIP-locked','Fixation-locked','Stimulus-locked'};
if nargin<4||isempty(xlim_ms), xlim_ms = [ERP.gip.times(1), ERP.gip.times(end)]; end

figure('Name',['ERPs — ' panel_title], 'Units','normalized', 'Position',[0.05 0.1 0.9 0.8]);
for r = 1:numel(locks)
    lk = locks{r}; t = ERP.(lk).times;
    for c = 1:numel(chan_labels)
        lab = chan_labels{c};
        ax = subplot(3, numel(chan_labels), (r-1)*numel(chan_labels)+c); hold(ax,'on'); grid(ax,'on');
        mS = ERP.(lk).std.(lab).mean; sS = ERP.(lk).std.(lab).sem;
        mD = ERP.(lk).dev.(lab).mean; sD = ERP.(lk).dev.(lab).sem;

        hLines = [];
        if ~isempty(mS)
            hS = plot(ax, t, mS, 'b-', 'LineWidth', 1.3); hLines(end+1)=hS;
            if ~isempty(sS), fill(ax, [t fliplr(t)], [mS-sS fliplr(mS+sS)], 'b','FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off'); end
        end
        if ~isempty(mD)
            hD = plot(ax, t, mD, 'r-', 'LineWidth', 1.3); hLines(end+1)=hD;
            if ~isempty(sD), fill(ax, [t fliplr(t)], [mD-sD fliplr(mD+sD)], 'r','FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off'); end
        end
        if ~isempty(hLines), lgd=legend(ax,hLines,{'Standard','Deviant'},'Location','northeast'); set(lgd,'Box','off'); end
        if r==1, title(ax, lab); end
        if c==1, ylabel(ax, sprintf('%s\n\\muV', row_titles{r})); end
        xline(ax,0,'k:','HandleVisibility','off'); xlim(ax,xlim_ms);
        if ~isempty(ylims), ylim(ax, ylims); end
        if r==numel(locks), xlabel(ax,'Time (ms)'); end
    end
end
sgtitle(panel_title,'FontWeight','bold');
end

function plot_erps_std_dev_multichan_with_stats(ERP, panel_title, chan_labels, xlim_ms, alpha_raw, q_fdr, mode, fdr_win)
% Adds gray (raw p<alpha) + black (FDR q<q_fdr) dots in fdr_win.
if nargin<4||isempty(xlim_ms), xlim_ms = [ERP.gip.times(1), ERP.gip.times(end)]; end
if nargin<5||isempty(alpha_raw), alpha_raw=0.05; end
if nargin<6||isempty(q_fdr), q_fdr=0.05; end
if nargin<7||isempty(mode), mode='trial_level'; end
if nargin<8||isempty(fdr_win), fdr_win=xlim_ms; end

ylims = [-4 4]; % e.g., [-4 4]
locks = {'gip','fix','stim'}; row_titles = {'GIP-locked','Fixation-locked','Stimulus-locked'};

figure('Name',['ERPs + FDR — ' panel_title], 'Units','normalized','Position',[0.05 0.1 0.9 0.8]);
for r = 1:numel(locks)
    lk = locks{r}; t = ERP.(lk).times; T = numel(t);
    for c = 1:numel(chan_labels)
        lab = chan_labels{c};
        ax = subplot(3, numel(chan_labels), (r-1)*numel(chan_labels)+c); hold(ax,'on'); grid(ax,'on');

        mS = ERP.(lk).std.(lab).mean; sS = ERP.(lk).std.(lab).sem;
        mD = ERP.(lk).dev.(lab).mean; sD = ERP.(lk).dev.(lab).sem;

        hLines = [];
        if ~isempty(mS)
            hS = plot(ax, t, mS, 'b-', 'LineWidth', 1.3); hLines(end+1)=hS;
            if ~isempty(sS), fill(ax,[t fliplr(t)],[mS-sS fliplr(mS+sS)],'b','FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off'); end
        end
        if ~isempty(mD)
            hD = plot(ax, t, mD, 'r-', 'LineWidth', 1.3); hLines(end+1)=hD;
            if ~isempty(sD), fill(ax,[t fliplr(t)],[mD-sD fliplr(mD+sD)],'r','FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off'); end
        end
        if ~isempty(hLines), lgd=legend(ax,hLines,{'Standard','Deviant'},'Location','northeast'); set(lgd,'Box','off'); end
        if r==1, title(ax, lab); end
        if c==1, ylabel(ax, sprintf('%s\n\\muV', row_titles{r})); end
        xline(ax,0,'k:','HandleVisibility','off'); xlim(ax,xlim_ms);
        if ~isempty(ylims), ylim(ax, ylims); end
        if r==numel(locks), xlabel(ax,'Time (ms)'); end

        % --- Stats (pointwise) ---
        Xs = ERP.trials.(lk).std.(lab); Xd = ERP.trials.(lk).dev.(lab);
        if isempty(Xs) || isempty(Xd), continue; end

        p_vals = nan(1,T);
        switch lower(mode)
            case 'trial_level'
                for tt = 1:T
                    xs = Xs(:,tt); xs = xs(~isnan(xs));
                    xd = Xd(:,tt); xd = xd(~isnan(xd));
                    if numel(xs)>=2 && numel(xd)>=2
                        [~, p_vals(tt)] = ttest2(xs, xd, 'Vartype','unequal');
                    end
                end
            otherwise
                error('Only ''trial_level'' is supported here (no file IDs in epochs).');
        end

        inwin = (t >= fdr_win(1) & t <= fdr_win(2));
        raw_mask = false(1,T); raw_mask(inwin) = p_vals(inwin) < alpha_raw;

        fdr_mask = false(1,T);
        p_win = p_vals(inwin); idx_local = find(inwin); good = ~isnan(p_win);
        if any(good)
            if exist('mafdr','file') == 2
                qvals = mafdr(p_win(good), 'BHFDR', true);
                sig = qvals < q_fdr;
            else
                % Fallback: use BH-adjusted p via fdr_bh; mark adj_p < q_fdr
                [~,~,adj_p] = fdr_bh(p_win(good), q_fdr, 'pdep', 'no');
                sig = adj_p < q_fdr;
            end
            % Place mask back into full-time vector
            full_idx = idx_local(good);
            fdr_mask(full_idx(sig)) = true;
        end

        yb = ylim(ax); ybase = yb(1);
        scatter(ax, t(raw_mask), repmat(ybase,1,sum(raw_mask)), 8, [0.7 0.7 0.7], 'filled', 'HandleVisibility','off');
        scatter(ax, t(fdr_mask), repmat(ybase,1,sum(fdr_mask)), 10, 'k', 'filled', 'HandleVisibility','off');
    end
end
sgtitle(sprintf('%s  (raw p<%.2f, FDR q<%.2f in %d–%d ms)', panel_title, alpha_raw, q_fdr, fdr_win(1), fdr_win(2)), 'FontWeight','bold');
end

function [h, crit_p, adj_p, sorted_p] = fdr_bh(pvals,q,method,report)
% Benjamini-Hochberg FDR, fallback if mafdr unavailable
if nargin < 2 || isempty(q), q = 0.05; end
if nargin < 3 || isempty(method), method = 'pdep'; end
if nargin < 4, report = 'no'; end
p = pvals(:); n = numel(p);
[sorted_p, sort_ids] = sort(p);
switch lower(method)
    case 'pdep', thresh = (1:n)'/n * q;
    case 'dep',  thresh = (1:n)'/n * q / sum(1./(1:n));
    otherwise, error('method must be ''pdep'' or ''dep''');
end
wtd_p = n * sorted_p ./ (1:n)'; adj_temp = cummin(wtd_p(end:-1:1)); adj_temp = adj_temp(end:-1:1);
adj_p = zeros(n,1); adj_p(sort_ids) = min(1, adj_temp);
h = sorted_p <= thresh; crit_p = any(h) * max(sorted_p(h));
if strcmpi(report,'yes')
    if any(h), fprintf('FDR threshold: %.4f\n', crit_p);
    else, fprintf('No significant results at q=%.2f\n', q); end
end
end

function debug_print_counts(src_name, side, labels, tAfter, rule, apply_flag)
    % labels: 1 x trials cell array ('standard'/'deviant')
    % tAfter: 1 x trials (ms)
    labels = labels(:);
    is_std0 = strcmpi(labels,'standard');
    is_dev0 = strcmpi(labels,'deviant');

    fprintf('\n[%s | %s]\n', src_name, side);
    fprintf('  Raw trials    : std=%d, dev=%d, total=%d\n', sum(is_std0), sum(is_dev0), numel(labels));

    if apply_flag
        keep = true(size(tAfter(:)));
        if isfield(rule,'min_ms') && ~isempty(rule.min_ms), keep = keep & (tAfter(:) >= rule.min_ms); end
        if isfield(rule,'max_ms') && ~isempty(rule.max_ms), keep = keep & (tAfter(:) <= rule.max_ms); end

        is_std = is_std0 & keep;
        is_dev = is_dev0 & keep;
        fprintf('  After trimming: std=%d, dev=%d (min=%s, max=%s)\n', ...
            sum(is_std), sum(is_dev), mat2str(rule.min_ms), mat2str(rule.max_ms));
    else
        fprintf('  No trimming applied to this lock.\n');
    end
end

