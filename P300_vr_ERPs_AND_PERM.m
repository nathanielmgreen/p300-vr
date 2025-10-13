clc; clear; close all;

%% =========================
%  User params (edit here)
% ==========================
mat_file   = 'epochs_data_08302025_new.mat';   % contains "epochs"
channels   = {'Fz','Pz','Oz','O1','O2'};       % channel labels to analyze
baseline_ms= [-300 0];                          % baseline for ERP plotting
xlim_ms    = [-300 1000];                      % plot x-limit
win_ms     = [0 400];                          % test window (ms)
stats_mode = 'trial_level';                    % 'trial_level' | 'paired_by_file' (paired_by_file not used unless you add IDs)
tail       = 'both';                           % 'both' | 'right' | 'left'
n_perm     = 5000;                             % permutations
use_cluster= false;                            % false: max-T; true: cluster-based
cluster_forming_p = 0.05;                      % only if use_cluster=true
rng(17);                                       % reproducible

% If you want to ignore your trimming policies for this dataset, just don’t
% filter trials here (this script uses whatever is in epochs.*.inner/outer_*).

%% =========================
%  Load
% ==========================
S = load(mat_file);
assert(isfield(S,'epochs'), 'File must contain struct "epochs".');
E = S.epochs;

% Locks we’ll analyze
locks = {'stimulus_onset','gip_onset','fix_onset'};
lock_labels = {'Stimulus-locked','GIP-locked','Fixation-locked'};

% Convenience: build channel index map
chan_labels = {E.chanlocs.labels};
chan_idx = cellfun(@(nm) find(strcmpi(chan_labels, nm),1), channels, 'UniformOutput', false);
if any(cellfun(@isempty, chan_idx))
    missing = channels(cellfun(@isempty, chan_idx));
    error('Channels not found: %s', strjoin(missing, ', '));
end
chan_idx = cell2mat(chan_idx);

%% =========================
%  Run for INNER and OUTER
% ==========================
sides = {'inner','outer'};
side_labels = struct('inner','Eye Shifting (Inner)','outer','Head Turning (Outer)');

for s = 1:numel(sides)
    side = sides{s};
    figure('Name',['Permutation ERPs — ' side_labels.(side)], ...
           'Units','normalized','Position',[0.05 0.06 0.9 0.84],'Color','w');

    % 3 rows (locks) × N channels columns
    for li = 1:numel(locks)
        L = locks{li};
        % Pull data & meta for this lock/side
        block = E.(L);
        data   = block.([side '_data']);       % trials × channels × time  OR  channels × time × trials  (we'll detect)
        stimuli= block.([side '_stimuli']);    % 1×trials cellstr (e.g., 'standard'/'deviant')
        tvec   = block.([side '_timeAfterStim']); % 1×time (ms)
        
        if isstring(stimuli) || ischar(stimuli)
            stimuli = cellstr(stimuli);
        end
        if iscategorical(stimuli)
            stimuli = cellstr(stimuli);
        end
        stimuli = stimuli(:).';
        nTrials = numel(stimuli);
        
        % trial masks per condition
        is_std = strcmpi(stimuli,'standard');
        is_dev = strcmpi(stimuli,'deviant');

        % For each channel, extract trials×time, do baseline (for plot), run permutation in win
        for ci = 1:numel(chan_idx)
            ch_i = chan_idx(ci);

            % Extract matrices: trials × time (handles common 3D layouts)
            X = trials_by_time(data, ch_i);
            
            X_std = X(is_std, :);
            X_dev = X(is_dev, :);

            % --- permutation stats (no baseline needed for testing if already mean-centered; using raw is fine)
            idx_win = tvec >= win_ms(1) & tvec <= win_ms(2);
            if ~any(idx_win)
                warning('%s/%s: test window not in time vector. Skipping.', side, L);
                raw_mask=false(size(tvec)); corr_mask=raw_mask;
            else
                switch lower(stats_mode)
                    case 'trial_level'
                        if use_cluster
                            corr_mask_win = cluster_perm_time(X_std(:,idx_win), X_dev(:,idx_win), @tt2_vec, n_perm, tail, cluster_forming_p);
                        else
                            corr_mask_win = perm_independent_maxT(X_std(:,idx_win), X_dev(:,idx_win), n_perm, tail);
                        end
                        corr_mask = false(size(tvec)); corr_mask(idx_win) = corr_mask_win;
                        % raw p<.05 (for gray dots)
                        p_raw = nan(1,numel(tvec));
                        for tt = find(idx_win)
                            [~,p_raw(tt)] = ttest2(X_std(:,tt), X_dev(:,tt),'Vartype','unequal');
                        end
                        raw_mask = p_raw < 0.05;

                    case 'paired_by_file'
                        error('paired_by_file requires subject/file IDs; not present in this epochs struct.');
                    otherwise
                        error('stats_mode must be trial_level or paired_by_file');
                end
            end

            % --- plotting (baseline for visualization)
            Xs = baseline_trials(X_std, tvec, baseline_ms);
            Xd = baseline_trials(X_dev, tvec, baseline_ms);
            mS = mean(Xs,1); sS = std(Xs,0,1)/sqrt(size(Xs,1));
            mD = mean(Xd,1); sD = std(Xd,0,1)/sqrt(size(Xd,1));

            % subplot index
            ax = subplot(numel(locks), numel(channels), (li-1)*numel(channels)+ci);
            hold(ax,'on'); grid(ax,'on');

            plot(tvec, mS, 'b','LineWidth',1.4);
            plot(tvec, mD, 'r','LineWidth',1.4);
            fill([tvec fliplr(tvec)], [mS-sS fliplr(mS+sS)], 'b','FaceAlpha',0.12,'EdgeColor','none');
            fill([tvec fliplr(tvec)], [mD-sD fliplr(mD+sD)], 'r','FaceAlpha',0.12,'EdgeColor','none');

            yl = ylim(ax); base = yl(1);
            % raw (uncorrected) gray
            scatter(tvec(raw_mask), repmat(base,1,sum(raw_mask)), 8, [0.7 0.7 0.7],'filled');
            % permutation-corrected black
            scatter(tvec(corr_mask), repmat(base,1,sum(corr_mask)), 10, 'k','filled');

            xline(0,'k:'); xlim(xlim_ms);
            if ci==1, ylabel([lock_labels{li} '  (\muV)']); end
            title(sprintf('%s — %s', channels{ci}, lock_labels{li}));
            if li==numel(locks)
                xlabel('Time (ms)');
            end
        end
    end

    % Single legend outside grid (create dummy handles)
    han = axes('Position',[0 0 1 1],'Visible','off');
    legend(han, {'Standard','Deviant','Raw p<.05','Perm FWER'}, 'Location','southoutside','Orientation','horizontal');
    sgtitle(['Permutation ERPs — ' side_labels.(side)]);
end

%% =========================
%        HELPERS
% =========================
function X = trials_by_time(data3, ch_i)
% Returns trials × time for channel ch_i from a 3D array.
% Supports: (trials × channels × time) OR (channels × time × trials) OR (time × channels × trials)
    sz = size(data3);
    if numel(sz) ~= 3
        error('Data must be 3D, got size %s', mat2str(sz));
    end
    % try common layouts
    if sz(2) >= ch_i && sz(3) >= 2 && sz(1) > 2
        % trials × channels × time
        X = squeeze(data3(:, ch_i, :));
    elseif sz(1) >= ch_i && sz(2) >= 2 && sz(3) > 2
        % channels × time × trials
        X = squeeze(data3(ch_i, :, :))';
    elseif sz(2) >= ch_i && sz(1) >= 2 && sz(3) > 2
        % time × channels × trials
        X = squeeze(data3(:, ch_i, :))';
    else
        error('Could not resolve data layout for channel index %d', ch_i);
    end
    % ensure trials × time
    if size(X,1) < size(X,2) && size(X,1) < 5
        % heuristic: if very few rows, probably transposed — keep as is
    end
end

function Xb = baseline_trials(X, tvec, base_ms)
% subtract per-trial mean over baseline window
    mask = tvec>=base_ms(1) & tvec<=base_ms(2);
    if ~any(mask) || isempty(X)
        Xb = X; return;
    end
    bl = mean(X(:,mask),2);
    Xb = X - bl;
end

%% ---------- Permutation engines ----------
function sig_mask = perm_independent_maxT(A, B, n_perm, tail)
% A,B: trials × time (independent groups). FWER via max-|T|.
    [na,T] = size(A); nb = size(B,1);
    % observed t
    t_obs = zeros(1,T);
    for t=1:T, t_obs(t)=welch_t(A(:,t),B(:,t)); end
    % maxT null
    pool = [A;B]; n = size(pool,1);
    maxT = zeros(n_perm,1);
    for p=1:n_perm
        idx = randperm(n);
        Ap = pool(idx(1:na),:);
        Bp = pool(idx(na+1:end),:);
        tp = zeros(1,T);
        for t=1:T, tp(t)=welch_t(Ap(:,t),Bp(:,t)); end
        maxT(p) = max_by_tail(tp, tail);
    end
    thr = prctile(maxT,95);
    sig_mask = pass_tail(t_obs, thr, tail);
end

function t = welch_t(x,y)
    nx=length(x); ny=length(y);
    vx=var(x,0);  vy=var(y,0);
    mx=mean(x);   my=mean(y);
    se = sqrt(vx/nx + vy/ny);
    if se==0, t=0; else, t=(mx-my)/se; end
end

function v = max_by_tail(tvec, tail)
    switch lower(tail)
        case 'both', v = max(abs(tvec));
        case 'right',v = max(tvec);
        case 'left', v = max(-tvec);
        otherwise, error('tail');
    end
end

function mask = pass_tail(t_obs, thr, tail)
    switch lower(tail)
        case 'both', mask = abs(t_obs) >= thr;
        case 'right',mask = t_obs >= thr;
        case 'left', mask = -t_obs >= thr;
        otherwise, error('tail');
    end
end

% Optional cluster-based (time) — uses cluster-forming p and max cluster mass.
function sig_mask = cluster_perm_time(A,B, tfun, n_perm, tail, p_form)
% A,B: trials × time; tfun: @tt2_vec (indep) or similar
    t_obs = tfun(A,B);
    df_approx = size(A,1)+size(B,1)-2;  % ok for Welch in large n
    p_obs = t_to_p(t_obs, df_approx, tail);
    clmask_obs = p_obs < p_form;
    masses_obs = cluster_masses(t_obs, clmask_obs);
    max_mass = zeros(n_perm,1);
    % permutations
    pool = [A;B]; nA = size(A,1); N=size(pool,1);
    for p=1:n_perm
        idx = randperm(N);
        Ap = pool(idx(1:nA),:);
        Bp = pool(idx(nA+1:end),:);
        tp = tfun(Ap,Bp);
        pp = t_to_p(tp, df_approx, tail);
        cl = pp < p_form;
        ms = cluster_masses(tp, cl);
        max_mass(p) = ifelse(isempty(ms),0,max(ms));
    end
    thr = prctile(max_mass,95);
    keep = masses_obs >= thr;
    sig_mask = false(1,length(t_obs));
    cl_idxs = clusters_indices(clmask_obs);
    for k=1:numel(cl_idxs)
        if keep(k), sig_mask(cl_idxs{k}) = true; end
    end
end

function tvec = tt2_vec(A,B)
    T=size(A,2); tvec=zeros(1,T);
    for t=1:T, tvec(t)=welch_t(A(:,t),B(:,t)); end
end

function p = t_to_p(t, df, tail)
    switch lower(tail)
        case 'both', p = 2*(1-tcdf(abs(t), df));
        case 'right',p = 1 - tcdf(t, df);
        case 'left', p = tcdf(t, df);
        otherwise, error('tail');
    end
end

function masses = cluster_masses(t_vec, bin_mask)
% sum |t| within each contiguous cluster
    idxs = clusters_indices(bin_mask);
    masses = zeros(1,numel(idxs));
    for k=1:numel(idxs), masses(k)=sum(abs(t_vec(idxs{k}))); end
end

function idxs = clusters_indices(mask)
    idxs = {};
    if ~any(mask), return; end
    d = diff([false mask false]);
    starts = find(d==1); ends = find(d==-1)-1;
    for i=1:numel(starts), idxs{end+1}=starts(i):ends(i); end
end

function y = ifelse(c,a,b), if c, y=a; else, y=b; end; end

function results = permtest_window_maxT(X_std, X_dev, tvec, win_ms, varargin)
% Permutation test in a window (e.g., 0–400 ms) following Groppe et al. (2011).
% Two-sample (independent) label-shuffle version; use sign-flip for paired.
% RETURNS: logical masks for corrected (FWER) and raw p<.05 across full tvec.
%
% Usage:
%   R = permtest_window_maxT(X_std, X_dev, tvec, [0 400], ...
%        'n_perm',5000,'tail','both','cluster',false,'p_form',0.05);

p = inputParser;
p.addParameter('n_perm', 5000);
p.addParameter('tail', 'both');       % 'both'|'right'|'left'
p.addParameter('cluster', false);     % false: max-|T|; true: cluster-mass
p.addParameter('p_form', 0.05);       % cluster-forming p (if cluster=true)
p.addParameter('alpha_raw', 0.05);    % for gray uncorrected dots
p.parse(varargin{:});
n_perm = p.Results.n_perm;
tail   = p.Results.tail;
use_cl = p.Results.cluster;
p_form = p.Results.p_form;
alpha_raw = p.Results.alpha_raw;

T  = numel(tvec);
idx_win = (tvec>=win_ms(1) & tvec<=win_ms(2));
if ~any(idx_win)
    error('permtest_window_maxT: window has no samples in tvec');
end

% ---------- Observed t across window ----------
t_obs = zeros(1,sum(idx_win));
ti = 0;
for t = find(idx_win)
    ti=ti+1;
    t_obs(ti) = welch_t(X_std(:,t), X_dev(:,t));  % two-sample Welch t
end

% ---------- Permutation null ----------
% Pool & shuffle labels under H0 (independent groups)
pool = [X_std; X_dev];
na   = size(X_std,1); N = size(pool,1);

if ~use_cl
    % ---- max-|T| (single-voxel) ----
    maxT = zeros(n_perm,1);
    for k=1:n_perm
        rp = randperm(N);
        Ap = pool(rp(1:na), idx_win);
        Bp = pool(rp(na+1:end), idx_win);
        tp = arrayfun(@(c) welch_t(Ap(:,c), Bp(:,c)), 1:size(Ap,2));
        maxT(k) = max_by_tail(tp, tail);
    end
    thr = prctile(maxT, 95);    % ~α = .05
    corr_mask_win = pass_tail(t_obs, thr, tail);
else
    % ---- cluster-based (time) with cluster-mass ----
    df_approx = size(X_std,1)+size(X_dev,1)-2;
    % observed clusters
    p_obs = t_to_p(t_obs, df_approx, tail);
    cf_mask_obs = p_obs < p_form;
    masses_obs  = cluster_masses(t_obs, cf_mask_obs);
    % permutation null = max cluster mass
    max_mass = zeros(n_perm,1);
    for k=1:n_perm
        rp = randperm(N);
        Ap = pool(rp(1:na), idx_win);
        Bp = pool(rp(na+1:end), idx_win);
        tp = arrayfun(@(c) welch_t(Ap(:,c), Bp(:,c)), 1:size(Ap,2));
        pp = t_to_p(tp, df_approx, tail);
        cf = pp < p_form;
        ms = cluster_masses(tp, cf);
        max_mass(k) = iff(isempty(ms), 0, max(ms));
    end
    thr = prctile(max_mass, 95);
    keep = masses_obs >= thr;
    corr_mask_win = false(1,sum(idx_win));
    clusters = clusters_indices(cf_mask_obs);
    for i=1:numel(clusters), if keep(i), corr_mask_win(clusters{i}) = true; end, end
end

% ---------- Expand masks to full time & also compute raw ----------
corr_mask = false(1,T); corr_mask(idx_win) = corr_mask_win;

p_raw = nan(1,T);
for t = find(idx_win)
    [~,p_raw(t)] = ttest2(X_std(:,t), X_dev(:,t), 'Vartype','unequal');
end
raw_mask = (p_raw < alpha_raw);

results.corr_mask = corr_mask;
results.raw_mask  = raw_mask;
results.t_obs     = t_obs;
results.thresh    = thr;
results.idx_win   = idx_win;
results.tail      = tail;
end

% ===== helpers =====
function t = welch_t(x,y)
nx=length(x); ny=length(y);
vx=var(x,0); vy=var(y,0);
mx=mean(x); my=mean(y);
se = sqrt(vx/nx + vy/ny);
if se==0, t=0; else, t=(mx-my)/se; end
end

function v = max_by_tail(tvec, tail)
switch lower(tail)
    case 'both', v = max(abs(tvec));
    case 'right',v = max(tvec);
    case 'left', v = max(-tvec);
end
end

function mask = pass_tail(t_obs, thr, tail)
switch lower(tail)
    case 'both', mask = abs(t_obs) >= thr;
    case 'right',mask = t_obs >= thr;
    case 'left', mask = -t_obs >= thr;
end
end

function p = t_to_p(t, df, tail)
switch lower(tail)
    case 'both', p = 2*(1-tcdf(abs(t), df));
    case 'right',p = 1 - tcdf(t, df);
    case 'left', p = tcdf(t, df);
end
end

function masses = cluster_masses(t_vec, bin_mask)
idxs = clusters_indices(bin_mask);
masses = zeros(1,numel(idxs));
for k=1:numel(idxs), masses(k)=sum(abs(t_vec(idxs{k}))); end
end

function idxs = clusters_indices(mask)
idxs = {};
if ~any(mask), return; end
d = diff([false mask false]);
s = find(d==1); e = find(d==-1)-1;
for i=1:numel(s), idxs{end+1} = s(i):e(i); end
end

function y = iff(c,a,b), if c, y=a; else, y=b; end, end
