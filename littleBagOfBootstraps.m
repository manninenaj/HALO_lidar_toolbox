function [lbob] = littleBagOfBootstraps(lbob, X, Y, Y_errors, method, display_msgs)
%littleBagofBootstraps function implements a procedure that uses
% bootstrapping together with subsampling, which enables a robust and
% computationally efficient assesment of an estimator and its quality.
% Kleiner et al., 2014. doi:10.1111/rssb.12050
% Adapted from Domino Data Science Platform project:
% - description: https://goo.gl/qMmZe1
% - examples: https://goo.gl/ZRxBrT
% - code: https://goo.gl/WsdDgM
%
% Usage:
% lbob = littleBagOfBootstraps(lbob, X, Y)
% lbob = littleBagOfBootstraps(lbob, X, Y, Y_errors)
%
% Inputs:
% - lbob                  struct
%       .score_func       function handle, a function used to score sampled 
%                         values. Special use: if string, can be either 
%                         'wstats', or 'wstats-weighted'
%       .agg_func         function handle, a function used to aggregate 
%                         separate estimates made in the inner loop.
%       .subsample_size   numeric, sample size in outer loop > 0 and <= 1.
%       .n_subsamples     numeric, number of sub-samples taken by the outer 
%                         loop. Three is minimum.
%       .n_trials         numeric, number of samples taken by the inner 
%                         loop. These samples are drawn repetitively from 
%                         a single sub-sample. Default is 30.
% - X                     numeric, vector or matrix. Can be empty, X = [];
% - Y                     numeric, vector or matrix
% - Y_errors              numeric, vector or matrix
%
% Outputs
% - lbob                  struct
%       .scores           estimates generated on each iteration of the 
%                         outer loop. These are used to estimate the upper 
%                         and lower bounds.
%       .best_estimate    final value, the best estimate calculated at the 
%                         outer loop
%       .standard_error   standard error of the estimate
%       .confidence_interval_low 
%                         confidence interval low, 25th percentile
%       .confidence_interval_high
%                         confidence interval high, 75th percentile
%
% Created 2018-02-01
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen@helsinki.fi

% Check inputs
if ~isstruct(lbob) || any(~isfield(lbob,'score_func') || ~isfield(lbob,'agg_func'))       
    error(sprintf(['The 1st input ''lbob'' needs to be a struct including, at least, '...
        'the fields ''score_func'' and ''agg_func''. \nType: help littleBagOfBootstraps']))   
end
if isa(lbob.score_func, 'function_handle') && ...
        ~(strcmp(func2str(lbob.score_func),'weightedMean') || ...
        strcmp(func2str(lbob.score_func),'weightedStddev') || ...
        strcmp(func2str(lbob.score_func),'weightedVariance') || ...
        strcmp(func2str(lbob.score_func),'weightedSkewness') || ...
        strcmp(func2str(lbob.score_func),'weightedKurtosis'))
    error(['If the field ''score_func'' in ''lbob'' is a function '...
        ' handle, it can be either a ''weightedMean'', ''weightedStddev'', '...
        ' ''weightedVariance'', ''weightedSkewness'', or ''weightedKurtosis''.'])
end
if (~isa(lbob.score_func, 'function_handle') && ...
        ~ischar(lbob.score_func)) || ...
        ~isa(lbob.score_func, 'function_handle') && ...
        ~any(strcmp(lbob.score_func,{'wstats','wstats-weighted'}))
    error(['If the field ''score_func'' in ''lbob'' is not a function '...
        ' handle, it can be either a ''wstats'' or a ''wstats-weighted'' string.'])
end
if isfield(lbob,'subsample_size')
    if ~isnumeric(lbob.subsample_size) || ...
        (lbob.subsample_size > 1 || lbob.subsample_size <= 0)
    error('The field ''subsample_size'' in ''lbob'' must be a number that is > 0 and <= 1.')
    end
else
    % Sample size in outer loop, b, n**0.67 used in literature.
    % If an integer it is taken as the number of rows to sample;
    % if less than 1, the sample size is taken as n**subsample_size.
    lbob.subsample_size = .67;
end
if isfield(lbob,'n_subsamples') 
    if ~isnumeric(lbob.n_subsamples) || (lbob.n_subsamples < 3)
        error('The field ''n_subsamples'' in ''lbob'' must be a number that is > 3.')
    end
else
    % The number of sub-samples taken by the outer loop. Three is the 
    % absolute minimum for simple distributions.
    lbob.n_subsamples = 3;
end
if isfield(lbob,'n_trials')
    if ~isnumeric(lbob.n_trials) || lbob.n_trials < 1
        error('The field ''n_trials'' in ''lbob'' must be a number that is >= 1.')
    end
else
    % The number of samples taken by the inner loop. These samples are 
    % drawn repetitively from a single sub-sample. Default is 30.
    lbob.n_trials = 30;
end
if ~isnumeric(Y(:)) || (~isvector(Y) && ~ismatrix(Y) && isscalar(Y))
    error('The input ''Y'' must be a numerical finite vector or matrix.')
end
if ~isempty(Y_errors) && (~isnumeric(Y_errors(:)) || (~isvector(Y_errors) && ...
        ~ismatrix(Y_errors) && isscalar(Y_errors)) || ...
        size(Y_errors,1)~=size(Y,1) || size(Y_errors,2)~=size(Y,2))
    error(['The input ''Y_errors'' must be a numerical finite vector or matrix,'...
        ' and has to have the same dimensions as ''Y''.'])
end

if nargin < 5
    method = 'weighted';
    display_msgs = true;
elseif nargin == 5
    if ~ischar(method) | strcmp(method,{'weighted','unweighted'})
        error('Input ''method'' can only be a string ''weighted'' or ''unweighted''.')
    end
    display_msgs = true;
end
if nargin < 6
    display_msgs = true;
elseif nargin == 6
    if ~islogical(display_msgs) & isscalar(display_msgs)
        error('The 6th input has to be scalar and logical true or false.')
    end
else
    error('Too many inputs.')
end

if ~isempty(X)
    lbob.x = X;
end
if ~isempty(Y)
    lbob.y = Y;
end
if isempty(lbob.x)
    lbob.x = zeros(size(lbob.y,1),1);
end

% Whether the scoring function uses value frequencies.
% This version uses measurement errors, which is not implemented by the
% original version --> dont use frequencies!
lbob.use_freqs = false;

lbob.y_errors = Y_errors;

Yunused = lbob.y;
Xunused = lbob.x;
Yunused_errors = lbob.y_errors;

subsample_size = lbob.subsample_size;
if subsample_size < 1
    subsample_size = round(size(lbob.y,1).^subsample_size);
end
n_trials = lbob.n_trials;
sample_size = size(lbob.y,1);
use_freqs = lbob.use_freqs;

if isa(lbob.score_func, 'function_handle')
    lbob.scores = zeros(lbob.n_subsamples, size(Y,2));
elseif any(strcmp(lbob.score_func,{'wstats','wstats-weighted'}))
    lbob.mean_scores = zeros(lbob.n_subsamples, size(Y,2));
    lbob.std_scores = zeros(lbob.n_subsamples, size(Y,2));
    lbob.var_scores = zeros(lbob.n_subsamples, size(Y,2));
%     lbob.scores_skewness = zeros(lbob.n_subsamples, size(Y,2));
%     lbob.scores_kurtosis = zeros(lbob.n_subsamples, size(Y,2));
else
    lbob = [];
    return
end

% Call inner bootstrap (little_boot) with different data n_subsample times.
for subset = 1:lbob.n_subsamples
    if size(Yunused,1) <= subsample_size
        Xunused = lbob.x;
        Yunused = lbob.y;
        Yunused_errors = lbob.y_errors;
    end
    if display_msgs
        fprintf('\n    LBOB: subset %s/%s ',num2str(subset),num2str(lbob.n_subsamples))
    end
    test_size = (size(Yunused,1) - subsample_size);
    [Xunused, Yunused, Yunused_errors] = train_test_split(Xunused,Yunused,Yunused_errors,test_size);

    if isa(lbob.score_func, 'function_handle')
        [lbob.scores(subset,:)] = lbob_little_boot(lbob,Xunused,Yunused,Yunused_errors,n_trials,sample_size,use_freqs,method,display_msgs);
    else
        [~,lbob.mean_scores(subset,:),lbob.std_scores(subset,:),lbob.var_scores(subset,:)] = ...
            lbob_little_boot(lbob,Xunused,Yunused,Yunused_errors,n_trials,sample_size,use_freqs,method,display_msgs);
%         scores_skewness(subset,:),scores_kurtosis(subset,:)
    end
end
if isa(lbob.score_func, 'function_handle')
%     lbob.scores = scores;
    lbob.best_estimate = nanmean(lbob.scores); % best estimate
    lbob.standard_error = nanstd(lbob.scores); % standard error
    lbob.confidence_interval_low = prctile(lbob.scores, 5);
    lbob.confidence_interval_high = prctile(lbob.scores, 95);
else
%     lbob.mean_scores = scores_mean;
%     lbob.std_scores = scores_std;
%     lbob.var_scores = scores_var;
%     lbob.skewness_scores = scores_skewness;
%     lbob.kurtosis_scores = scores_kurtosis;
    
    lbob.mean_best_estimate = nanmean(lbob.mean_scores); % best estimate
    lbob.std_best_estimate = nanmean(lbob.std_scores); % best estimate
    lbob.var_best_estimate = nanmean(lbob.var_scores); % best estimate
%     lbob.skewness_best_estimate = nanmean(scores_skewness); % best estimate
%     lbob.kurtosis_best_estimate = nanmean(scores_kurtosis); % best estimate
    
    lbob.mean_standard_error = nanstd(lbob.mean_scores); % standard error
    lbob.std_standard_error = nanstd(lbob.std_scores); % standard error
    lbob.var_standard_error = nanstd(lbob.var_scores); % standard error
%     lbob.skewness_standard_error = nanstd(scores_skewness); % standard error
%     lbob.kurtosis_standard_error = nanstd(scores_kurtosis); % standard error
    
    lbob.mean_confidence_interval_low = prctile(lbob.mean_scores, 5);
    lbob.std_confidence_interval_low = prctile(lbob.std_scores, 5);
    lbob.var_confidence_interval_low = prctile(lbob.var_scores, 5);
%     lbob.skewness_confidence_interval_low = prctile(scores_skewness, 5);
%     lbob.kurtosis_confidence_interval_low = prctile(scores_kurtosis, 5);
    
    lbob.mean_confidence_interval_high = prctile(lbob.mean_scores, 95);
    lbob.std_confidence_interval_high = prctile(lbob.std_scores, 95);
    lbob.var_confidence_interval_high = prctile(lbob.var_scores, 95);
%     lbob.skewness_confidence_interval_high = prctile(scores_skewness, 95);
%     lbob.kurtosis_confidence_interval_high = prctile(scores_kurtosis, 95);
end
% fprintf('\n')
end

function [scores,scores_mean,scores_std,scores_var] = ...
    lbob_little_boot(lbob,X,Y,Y_errors,n_trials,sample_size,use_freqs,method,display_msgs)
% subsample_score_skewness,subsample_score_kurtosis
% Implements the inner loop of little-bag-of-bootstraps.
%
%     Carries out n_trials of resampling with replacement and calculates
%     a score on each sample using the scoring function. Applies the
%     aggregation function to the scores, usually a mean, and returns
%     a single value.
%
% X: Feature matrix.
% Y: Vector of values or dependent variables.
% n_trials: Number of times to sample.
% sample_size: Effective sample size of each sample.
% use_freqs: If True then frequencies are used.

% freq = [];
if isa(lbob.score_func, 'function_handle')
    scores = zeros(n_trials, size(Y,2));
else
    scores_mean = zeros(n_trials, size(Y,2));
    scores_std = zeros(n_trials, size(Y,2));
    scores_var = zeros(n_trials, size(Y,2));
%     scores_skewness = zeros(n_trials, size(Y,2));
%     scores_kurtosis = zeros(n_trials, size(Y,2));
end
for t = 1:n_trials

%     if use_freqs % use frequencies
% %         subsample_score = []; % TBD!?!?!? or not...
%         return
% %         [~, Y_eval, Y_errors_eval, ~] = lbob_get_multinomial_sample_with_freqs(X, Y, Y_errors, sample_size);
%     else % do not use frequencies
        [~, Y_eval, Y_errors_eval] = lbob_get_multinomial_sample(X, Y, Y_errors, sample_size);
%     end
%     if use_freqs % use frequencies
% %         subsample_score = []; % TBD!?!?!? or not...
%         return
%     else % do not use frequencies
        if isa(lbob.score_func, 'function_handle')
                switch method
                    case 'unweighted'
                        scores(t,:) = lbob.score_func(Y_eval,Y_errors_eval,method);
                    case 'weighted'
                        scores(t,:) = lbob.score_func(Y_eval,Y_errors_eval,method);
                end
        else
            switch lbob.score_func
                case 'wstats'
                    scores_mean(t,:) = weightedMean(Y_eval,Y_errors_eval,method);
                    scores_std(t,:) = weightedStddev(Y_eval,Y_errors_eval,method);
                    scores_var(t,:) = scores_std(t,:).^2;   
% % %                     scores_var(t,:) = weightedVariance(Y_eval,Y_errors_eval,method); DOES NOT WORK
%                     scores_skewness(t,:) = weightedSkewness(Y_eval,Y_errors_eval,'unweighted');
%                     scores_kurtosis(t,:) = weightedKurtosis(Y_eval,Y_errors_eval,'unweighted');
                case 'wstats-weighted'
                    scores_mean(t,:) = weightedMean(Y_eval,Y_errors_eval,method);
                    scores_std(t,:) = weightedStddev(Y_eval,Y_errors_eval,method);
                    scores_var(t,:) = scores_std(t,:).^2;
% % %                     scores_var(t,:) = weightedVariance(Y_eval,Y_errors_eval,method); DOES NOT WORK
%                     scores_skewness(t,:) = weightedSkewness(Y_eval,Y_errors_eval);
%                     scores_kurtosis(t,:) = weightedKurtosis(Y_eval,Y_errors_eval);
                otherwise
            end
        end
%     end
    if display_msgs        
        fprintf('%3.3g%% ',t/n_trials*100)
    end
end
% average (default) the aggregated scores
if isa(lbob.score_func, 'function_handle')
    scores = lbob.agg_func(scores);
else
    scores = [];
    scores_mean = lbob.agg_func(scores_mean);
    scores_std = lbob.agg_func(scores_std);
    scores_var = lbob.agg_func(scores_var);
%     scores_skewness = lbob.agg_func(scores_skewness);
%     scores_kurtosis = lbob.agg_func(scores_kurtosis);
end
end

function [X,Y,Y_errors] = lbob_get_multinomial_sample(X, Y, Y_errors, sample_size)
% Draw a multinomial sample from a vector or feature matrix, vector combination.
%
%     Samples with replacement. Does not return a set of frequencies.
%
% X: Feature matrix.
% Y: Dependent variables or values.
% sample_size: Number of rows in returned arrays. Defaults to size of input vector Y.

if isempty(sample_size)
    sample_size = size(Y,1);
end
sizeY = size(Y,1);
rows = randi([1 sizeY], [sample_size 1]);
Y = Y(rows,:);
if isempty(Y_errors)
    Y_errors = [];
else
    Y_errors = Y_errors(rows,:);
end
if ~isempty(X)
    X = X(rows, :);
end
clearvars rows
end

function [x,y,y_error] = train_test_split(x,y,y_error,test_size,random_state)
% [x_train,x,y_train,y,y_train_error,y_error] = train_test_split(x,y,y_error,test_size,random_state)

if nargin<4
    test_size = 1/3;
end

if nargin>4
    rng('default');
    rng(random_state);
end

n = size(x,1);
% n_test = floor(test_size*n);
% n_test = test_size;
% n_train = n - test_size;

indp = randperm(n);
% ind_train = indp(1:n_train);
% ind_test = indp(n_train+1:n);

x = x(indp(1:n-test_size),:);
% x = x(indp(n-test_size+1:n),:);
y = y(indp(1:n-test_size),:);
% y = y(indp(n-test_size+1:n),:);
if isempty(y_error)
    y_error = [];
%     y_error = [];
else
    y_error = y_error(indp(1:n-test_size),:);
%     y_error = y_error(indp(n-test_size+1:n),:);
end
end


%% TBD etc..

% function [Xsamp, Ysamp, Ysamp_errors, freq] = lbob_get_multinomial_sample_with_freqs(X, Y, Y_errors, sample_size)
% % Draw a multinomial sample from a vector or feature matrix, vector combination.
% %
% %     Returns the rows sampled and frequencies with which they were sampled.
% %     Unsampled rows are not returned.
% %
% % X: Feature matrix.
% % Y: Dependent variables or values.
% % sample_size: Number of rows in returned arrays. Defaults to size of input vector Y.
% 
% if isempty(sample_size)
%     sample_size = size(Y,1); 
% end
%     n_rows = size(Y,1);
%     freq_by_row = mnrnd(sample_size,repmat(1./n_rows,n_rows,1));
%     rows = find(freq_by_row);
%     freq = freq_by_row(rows);
%     Ysamp = Y(rows,:);
%     if isempty(Y_errors)
%         Ysamp_errors = [];
%     else
%         Ysamp_errors = Y_errors(rows,:);
%     end
% if ~isempty(X)
%     Xsamp = X(rows, :);
% end
% end

