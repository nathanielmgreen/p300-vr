clc; clear; close all;

%% =============================
%            SETUP
% ==============================
mat_file    = 'epochs_data_08302025_new.mat';   % <-- your .mat
chan_labels = {'Cz','Pz','Oz','O1','O2'};       % grid columns
baseline_ms = [-300 0];                         % baseline window (ms)
xlim_ms     = [-300 800];                       % plot window

% analysis windows
stats_win_inner = [0 400];    % eye shifting (inner)
stats_win_outer = [200 600];  % head turning (outer)

alpha_raw    = 0.05;          % raw (uncorrected) p threshold (gray)
alpha_holm   = 0.05;          % Holm–Bonferroni FWER level (black)
stats_mode   = 'trial_level'; % epochs struct lacks per-file IDs

% --- Permutation settings ---
n_perm            = 5000;     % increase as needed
tail              = 'left';   % 'both' | 'right' | 'left'  (e.g., Deviant < Standard if 'left')
prefilter_top_pct = 0.05;     % 0.05 = top 5% |t| within stats_win before Holm; 0 disables

% --- Trimming policy knobs (kept but disabled per your request) ---
rules_inner.gip = struct('min_ms', 50,  'max_ms', 400);
rules_inner.fix = struct('min_ms', 80,  'max_ms', 700);
rules_outer.gip = struct('min_ms', 50,  'max_ms', 600);
rules_outer.fix = struct('min_ms', 120, 'max_ms', 1200);
apply_rules = struct('gip', false, 'fix', false, 'stim', false); % keep OFF

%% =============================
%        LOAD EPOCHS STRUCT
% ==============================
S = load(mat_file);
if ~isfield(S,'epochs'), error('File must contain variable "epochs".'); end
epochs = S.epochs;
times  = epochs.times(:)';  % actual ERP timestamps

% map requested channel labels to indices
chIdxMap = map_channel_indices(epochs.chanlocs, chan_labels);

%% =============================
%   BUILD ERPs (INNER / OUTER)
% ==============================
ERP_INNER = build_erps_from_epochs(epochs, 'inner', times, chIdxMap, chan_labels, baseline_ms, rules_inner, apply_rules);
ERP_OUTER = build_erps_from_epochs(epochs, 'outer', times, chIdxMap, chan_labels, baseline_ms, rules_outer, apply_rules);

%% =============================
%          PLOTTING
% ==============================
% Plain ERPs (no stats)
plot_erps_std_dev_multichan(ERP_INNER, 'Eye Shifting (Inner)', chan_labels, xlim_ms);
plot_erps_std_dev_multichan(ERP_OUTER, 'Head Turning (Outer)', chan_labels, xlim_ms);

% ERPs with stats:
%  gray  = raw parametric p<alpha_raw
%  black = significant after Permutation → Holm–Bonferroni (FWER), optionally prefiltered to top-|t|%
plot_erps_permHOLM(ERP_INNER, 'Eye Shifting (Inner)', chan_labels, xlim_ms, ...
                   stats_win_inner, alpha_raw, alpha_holm, stats_mode, n_perm, tail, prefilter_top_pct);

plot_erps_permHOLM(ERP_OUTER, 'Head Turning (Outer)', chan_labels, xlim_ms, ...
                   stats_win_outer, alpha_raw, alpha_holm, stats_mode, n_perm, tail, prefilter_top_pct);

disp('Done.');

%% =============================
%            FUNCTIONS
% ==============================

function ERP = build_erps_from_epochs(epochs, side, times, chIdxMap, chan_labels, baseline_ms, rules, apply_rules)
% side: 'inner' or 'outer'
locks_src = {'gip_onset','fix_onset','stimulus_onset'};
locks_out = {'gip','fix','stim'};
ERP = struct();

for L = 1:numel(locks_src)
    src = locks_src{L};
    out = locks_out{L};

    data   = epochs.(src).([side '_data']);        % ch x time x trials
    labels = epochs.(src).([side '_stimuli']);     % 1 x trials cell
    tAfter = epochs.(src).([side '_timeAfterStim']); %#ok<NASGU> % not used when apply_rules=false

    % normalize labels
    if ischar(labels) || isstring(labels), labels = cellstr(labels); end
    if iscategorical(labels),              labels = cellstr(labels); end
    labels = labels(:).';

    is_std = strcmpi(labels,'standard');
    is_dev = strcmpi(labels,'deviant');

    % baseline (per-trial per-channel)
    if ~isempty(baseline_ms)
        data = baseline_epochs(data, times, baseline_ms);
    end

    ERP.(out).times = times;
    for c = 1:numel(chan_labels)
        lab = chan_labels{c};
        idx = chIdxMap.(lab);
        if isempty(idx)
            ERP.trials.(out).std.(lab) = [];
            ERP.trials.(out).dev.(lab) = [];
            ERP.(out).std.(lab) = summarize_trials_struct([]);
            ERP.(out).dev.(lab) = summarize_trials_struct([]);
            continue;
        end
        X = squeeze(data(idx, :, :))';  % trials x time
        ERP.trials.(out).std.(lab) = X(is_std, :);
        ERP.trials.(out).dev.(lab) = X(is_dev, :);
        ERP.(out).std.(lab) = summarize_trials_struct(ERP.trials.(out).std.(lab));
        ERP.(out).dev.(lab) = summarize_trials_struct(ERP.trials.(out).dev.(lab));
    end
end
end

function chIdxMap = map_channel_indices(chanlocs, wanted_labels)
chIdxMap = struct();
have = arrayfun(@(x) string(x.labels), chanlocs, 'UniformOutput', true);
for k = 1:numel(wanted_labels)
    lab = wanted_labels{k};
    idx = find(strcmpi(have, lab), 1);
    if isempty(idx)
        warning('Channel %s not found.', lab);
        chIdxMap.(lab) = [];
    else
        chIdxMap.(lab) = idx;
    end
end
end

function data_bc = baseline_epochs(data, times, baseline_ms)
data_bc = data;
t = times(:)'; bmin = max(baseline_ms(1), t(1)); bmax = min(baseline_ms(2), t(end));
if ~(bmin < bmax)
    warning('baseline_epochs: invalid baseline window; skipping.'); return;
end
bmask = (t>=bmin & t<=bmax);
if ~any(bmask)
    warning('baseline_epochs: no baseline samples; skipping.'); return;
end
for tr = 1:size(data,3)
    bl = mean(data(:, bmask, tr), 2, 'omitnan');
    data_bc(:, :, tr) = data(:, :, tr) - bl;
end
end

function S = summarize_trials_struct(X)
if isempty(X)
    S.mean = []; S.sem = []; S.n = 0;
else
    S.mean = mean(X,1,'omitnan');
    S.sem  = std(X,0,1,'omitnan')/sqrt(size(X,1));
    S.n    = size(X,1);
end
end

%% ---------- PLOTTING: plain ----------
function plot_erps_std_dev_multichan(ERP, panel_title, chan_labels, xlim_ms)
ylims = []; % e.g., [-4 4]
locks = {'gip','fix','stim'}; row_titles = {'GIP-locked','Fixation-locked','Stimulus-locked'};
if nargin<4||isempty(xlim_ms), xlim_ms = [ERP.gip.times(1), ERP.gip.times(end)]; end

figure('Name',['ERPs — ' panel_title], 'Units','normalized','Position',[0.05 0.1 0.9 0.8]);
for r=1:numel(locks)
    lk = locks{r}; t = ERP.(lk).times;
    for c=1:numel(chan_labels)
        lab = chan_labels{c};
        ax = subplot(3, numel(chan_labels), (r-1)*numel(chan_labels)+c); hold(ax,'on'); grid(ax,'on');
        mS = ERP.(lk).std.(lab).mean; sS = ERP.(lk).std.(lab).sem;
        mD = ERP.(lk).dev.(lab).mean; sD = ERP.(lk).dev.(lab).sem;

        if ~isempty(mS)
            plot(ax, t, mS, 'b-', 'LineWidth',1.3);
            if ~isempty(sS), fill(ax,[t fliplr(t)],[mS-sS fliplr(mS+sS)],'b','FaceAlpha',0.15,'EdgeColor','none'); end
        end
        if ~isempty(mD)
            plot(ax, t, mD, 'r-', 'LineWidth',1.3);
            if ~isempty(sD), fill(ax,[t fliplr(t)],[mD-sD fliplr(mD+sD)],'r','FaceAlpha',0.15,'EdgeColor','none'); end
        end
        if r==1, title(ax, lab); end
        if c==1, ylabel(ax, sprintf('%s\n\\muV', row_titles{r})); end
        xline(ax,0,'k:'); xlim(ax, xlim_ms);
        if ~isempty(ylims), ylim(ax, ylims); end
        if r==numel(locks), xlabel(ax,'Time (ms)'); end
    end
end
sgtitle(panel_title,'FontWeight','bold');
end

%% ---------- PLOTTING: Permutation → Holm–Bonferroni ----------
function plot_erps_permHOLM(ERP, panel_title, chan_labels, xlim_ms, stats_win, alpha_raw, alpha_holm, mode, n_perm, tail, prefilter_top_pct)
% ERPs with stats:
%  - gray dots/bands  = raw parametric p<alpha_raw in stats_win
%  - black dots/bands = significant after Permutation→Holm–Bonferroni (FWER), optionally prefiltered

locks = {'gip','fix','stim'}; row_titles = {'GIP-locked','Fixation-locked','Stimulus-locked'};
ylims = []; % e.g., [-4 4]
if nargin<4||isempty(xlim_ms), xlim_ms = [ERP.gip.times(1), ERP.gip.times(end)]; end

% Pillar style
pillar_width_ms   = [];    % [] → 4*median(diff(t)) per axis
alpha_gray_pillar = 0.18;  % opacity for raw pillars
alpha_blk_pillar  = 0.22;  % opacity for perm→Holm pillars
color_gray        = [0.7 0.7 0.7];
color_black       = [0   0   0  ];

figure('Name',['ERPs + Permutation→Holm (FWER) — ' panel_title], ...
       'Units','normalized','Position',[0.05 0.1 0.9 0.8]);

for r=1:numel(locks)
    lk = locks{r}; t = ERP.(lk).times; T=numel(t);
    inwin = (t>=stats_win(1) & t<=stats_win(2));
    for c=1:numel(chan_labels)
        lab = chan_labels{c};
        ax = subplot(3, numel(chan_labels), (r-1)*numel(chan_labels)+c); hold(ax,'on'); grid(ax,'on');

        Xs = ERP.trials.(lk).std.(lab); Xd = ERP.trials.(lk).dev.(lab);
        if isempty(Xs) || isempty(Xd)
            if r==1, title(ax, lab); end
            if c==1, ylabel(ax, sprintf('%s\n\\muV', row_titles{r})); end
            xline(ax,0,'k:'); xlim(ax, xlim_ms);
            continue;
        end

        % Means/SEMs
        mS = mean(Xs,1,'omitnan'); sS=std(Xs,0,1,'omitnan')/sqrt(size(Xs,1));
        mD = mean(Xd,1,'omitnan'); sD=std(Xd,0,1,'omitnan')/sqrt(size(Xd,1));

        % Plot waveforms first
        hS = plot(ax, t, mS, 'b', 'LineWidth',1.3, 'DisplayName','Standard');
        hD = plot(ax, t, mD, 'r', 'LineWidth',1.3, 'DisplayName','Deviant');
        if ~isempty(sS), fill(ax,[t fliplr(t)], [mS-sS fliplr(mS+sS)], 'b','FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off'); end
        if ~isempty(sD), fill(ax,[t fliplr(t)], [mD-sD fliplr(mD+sD)], 'r','FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off'); end

        % ---------- Raw parametric p (gray) ----------
        p_raw = nan(1,T);
        if strcmpi(mode,'trial_level')
            for tt=1:T
                xs = Xs(:,tt); xs=xs(~isnan(xs));
                xd = Xd(:,tt); xd=xd(~isnan(xd));
                if numel(xs)>=2 && numel(xd)>=2
                    [~,p_raw(tt)] = ttest2(xs, xd, 'Vartype','unequal');
                end
            end
        else
            error('Only trial_level supported with epochs struct.');
        end
        raw_mask = false(1,T); raw_mask(inwin) = p_raw(inwin) < alpha_raw;

        % ---------- Permute per-time, then Holm step-down within window (optionally prefiltered) ----------
        R = permHOLM_pointwise(Xs, Xd, t, stats_win, n_perm, tail, alpha_holm, prefilter_top_pct);

        % Lock axes before drawing pillars/dots
        xline(ax,0,'k:','HandleVisibility','off');
        xlim(ax, xlim_ms);
        if ~isempty(ylims), ylim(ax, ylims); end
        ax.XLimMode = 'manual'; ax.YLimMode = 'manual';
        yb = ylim(ax); ybase = yb(1);

        % dynamic pillar width
        if isempty(pillar_width_ms)
            if numel(t)>1, dt = median(diff(t)); else, dt = 1; end
            width_ms = 4*dt;
        else
            width_ms = pillar_width_ms;
        end

        % Translucent pillars
        draw_vertical_bands(ax, t(raw_mask), yb, color_gray,  alpha_gray_pillar, width_ms);
        draw_vertical_bands(ax, t(R.holm_mask), yb, color_black, alpha_blk_pillar, width_ms);

        % Dots at baseline
        scatter(ax, t(raw_mask),    repmat(ybase,1,sum(raw_mask)),    8,  color_gray, 'filled', 'HandleVisibility','off');
        scatter(ax, t(R.holm_mask), repmat(ybase,1,sum(R.holm_mask)), 10, 'k',         'filled', 'HandleVisibility','off');

        % Cosmetics
        if r==1, title(ax, lab); end
        if c==1, ylabel(ax, sprintf('%s\n\\muV', row_titles{r})); end
        if r==numel(locks), xlabel(ax,'Time (ms)'); end
        set(ax,'Layer','top');
        if r==1
            legend(ax, [hS hD], {'Standard','Deviant'}, 'Location','northeast','Box','off');
        end
    end
end

sgtitle(sprintf('%s  (raw p<%.2f, Perm→Holm FWER α=%.2f in %d–%d ms%s) — bands = translucent verticals', ...
    panel_title, alpha_raw, alpha_holm, stats_win(1), stats_win(2), ...
    tern(prefilter_top_pct>0, sprintf(', prefilter top %.0f%% |t|', 100*prefilter_top_pct), '')), ...
    'FontWeight','bold');
end

%% ---------- Permutation p’s → Holm–Bonferroni step-down ----------
function out = permHOLM_pointwise(Xa, Xb, tvec, win_ms, n_perm, tail, alpha_holm, top_pct)
% Per-timepoint permutation p-values within window, optional prefilter by top-|t|,
% then Holm–Bonferroni across tested points.

T  = numel(tvec);
win_idx = find(tvec>=win_ms(1) & tvec<=win_ms(2));
if isempty(win_idx), error('permHOLM_pointwise: empty window.'); end
na = size(Xa,1); nb = size(Xb,1); %#ok<NASGU>
pool = [Xa; Xb];
N   = size(pool,1);

% Observed t per time (Welch)
t_obs = nan(1,T);
for tt = win_idx(:)'
    t_obs(tt) = welch_t(Xa(:,tt), Xb(:,tt));
end

% Permuted null per time (label-shuffle), one-sided / two-sided p
p_perm = nan(1,T);
for tt = win_idx(:)'
    tperm = zeros(n_perm,1);
    for k=1:n_perm
        rp = randperm(N);
        Ap = pool(rp(1:na), tt);
        Bp = pool(rp(na+1:end), tt);
        tperm(k) = welch_t(Ap, Bp);
    end
    switch lower(tail)
        case 'both'
            p = mean(abs(tperm) >= abs(t_obs(tt)));
        case 'right'
            p = mean(tperm >= t_obs(tt));
        case 'left'
            p = mean(tperm <= t_obs(tt));
        otherwise
            error('tail must be both/right/left');
    end
    if p==0, p = 1/(n_perm+1); end   % small-sample continuity
    p_perm(tt) = p;
end

% Optional prefilter by top-|t| within window
use_idx = win_idx;
if top_pct > 0
    [~, order] = sort(abs(t_obs(win_idx)), 'descend');
    kkeep = max(1, round(top_pct * numel(win_idx)));
    use_idx = win_idx(order(1:kkeep));
end

% Holm–Bonferroni across the tested set
holm_mask = false(1,T);
if ~isempty(use_idx)
    p_sub = p_perm(use_idx);
    good  = ~isnan(p_sub);
    if any(good)
        pass = holm_stepdown(p_sub(good), alpha_holm);
        holm_mask(use_idx(good(pass))) = true;
    end
end

out.holm_mask = holm_mask;
end

%% ---------- Multiple-comparison helpers ----------
function pass = holm_stepdown(p, alpha)
% Holm–Bonferroni step-down procedure.
% Input p (vector), alpha scalar. Output pass (logical) same size as p.
p = p(:); m = numel(p);
[sp, idx] = sort(p, 'ascend');
pass_sorted = false(m,1);
% find largest k such that sp(j) <= alpha/(m-j+1) for all j=1..k
kstar = 0;
for j = 1:m
    if sp(j) <= alpha/(m - j + 1)
        kstar = j;
    else
        break
    end
end
if kstar > 0
    pass_sorted(1:kstar) = true; % reject H1..Hk*
end
pass = false(size(p));
pass(idx) = pass_sorted;
end

%% ---------- Small helpers ----------
function t = welch_t(x,y)
x = x(~isnan(x)); y = y(~isnan(y));
nx=numel(x); ny=numel(y);
if nx<2 || ny<2, t = 0; return; end
vx=var(x,0); vy=var(y,0);
mx=mean(x); my=mean(y);
se = sqrt(vx/nx + vy/ny);
if se==0, t=0; else, t=(mx-my)/se; end
end

function y = tern(c,a,b), if c, y=a; else, y=b; end, end

function draw_vertical_bands(ax, x_positions, ylims, rgb, alpha, width_ms)
% Draw semi-transparent full-height vertical bands centered at x_positions.
if isempty(x_positions), return; end
hb = gobjects(0);
for k = 1:numel(x_positions)
    x0 = x_positions(k) - width_ms/2;
    x1 = x_positions(k) + width_ms/2;
    hb(end+1) = patch(ax, [x0 x1 x1 x0], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
                      rgb, 'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility','off'); %#ok<AGROW>
end
try uistack(hb, 'bottom'); catch, end
end
