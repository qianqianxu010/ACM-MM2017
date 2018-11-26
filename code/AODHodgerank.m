function [score, output] = AODHodgerank(data, options)
%--------------------------------------------------------------------------
% AODHodgerank.m: Adaptive outlier detection with HodgeRank
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Adaptive outlier detection with HodgeRank
%
% USAGE: 
%    [score, num_outlier, outlier_detect] = AODHodgerank(data, options)
%    [score, num_outlier, outlier_detect] = AODHodgerank(data)
%    score = AODHodgerank(data)
%
% INPUT ARGUMENTS:
% data        Input matrix with paired comparison, the matrix has 
%             2 columns, and for each row, the rank of the first column
%             is higher than the second column for this comparison.
% options     The parameters are stored here, usually just the default 
%             parameters are good.
% options.intercept     The gap between the scores to be considered
%                       different.
% options.alpha         For the first step, estimate the number of outliers
%                       using the number of outliers based on HodgeRank.
%                       Default is 0.9
% options.beta1         Updating the estimation of the number of outliers.
%                       Default is 0.9, this number should be < 1.
% options.beta2         Updating the estimation of the number of outliers.
%                       Default is 1.01, this number should be > 1.
% OUTPUT ARGUMENTS:
% score                 The final scores which are integers from 0 to n-1 for adaHodgeRank 
%                       or others fro HodgeRank.
% output.
% num_outlier           Output of all store numbers of outliers: estimation, 
%                       based on the output score, two corrections.
% outlier_detect        The dection of outlier for each iteration.
% adjust                The adjust at the last step
% Z0                    The matrix of comparison
%
% LICENSE: GPL-2
%
% DATE: 2014-3-14
%
% AUTHORS: Ming Yan
%
% REFERENCES:
%
% SEE ALSO:
% Hodgerank, AdaHodgerank
%
%%
n               = max(data(:));
numsample       = size(data,1);
if nargin < 2
    options.intercept = 1;
end

if ~isfield(options,'alpha')
    options.alpha = 0.9;
end
if ~isfield(options,'beta1')
    options.beta1 = 0.9;
end
if ~isfield(options,'beta2')
    options.beta2 = 1.01;
end
if ~isfield(options,'adaHodge')
    options.adaHodge = false;
end
if ~isfield(options,'max_iter')
    options.max_iter = 30;
end


[score, Z0] = Hodgerank(data); % Use Hodegrank to obtain first first estimation

diff  = max(score(data(:,2))' - score(data(:,1))', 0);

[diff_sort index] = sort(diff,'descend');
%count=1;


num_outlier2    = sum(diff > 0);  % num_outlier2 is the outlier number for each output
num_outlier1    = floor(num_outlier2 * options.alpha); % estimation of the outlier number
num_outlier3    = sum(diff_sort > diff_sort(num_outlier1)); % correction of the outlier number
num_outlier4    = sum(diff_sort >= diff_sort(num_outlier1)); % correction 2 of the outlier number

if num_outlier1*(options.beta2 - 1) < 1
    options.beta2 = 1.5/num_outlier1 + 1;
    keyboard
elseif num_outlier1*(options.beta2 - 1) > 30
    options.beta2 = 25/num_outlier1 + 1;
end

num_outlier(1,:)    = [num_outlier1, num_outlier2, num_outlier3, num_outlier4];
num_outlier_B   = num_outlier2;
score_B         = score';
score           = ones(n,1);  % initial score 
for iter = 1:options.max_iter
    data2       = data(index(num_outlier1 + 1:numsample),:); % the data without the detected outliers
    outlier_detect(:,iter) = diff > diff_sort(num_outlier1);        
    score       = Hodgerank(data2, score)'; % Update score based on the "true" data
    
    diff         = max(score(data(:,2))' - score(data(:,1))', 0);

    [diff_sort index]   = sort(diff,'descend');
    
    num_outlier2    = sum(diff > 0);  % num_outlier2 is the outlier number for each output
    if num_outlier2 < num_outlier_B
        num_outlier_B = num_outlier2;
        score_B = score;
    end
    num_outlier3    = sum(diff_sort > diff_sort(num_outlier1)); % correction of the outlier number
    num_outlier4    = sum(diff_sort >= diff_sort(num_outlier1)); % correction 2 of the outlier number    
    if num_outlier2 < num_outlier(iter,2) % Better score, increase estimation
        num_outlier1    = max(floor(num_outlier2 * options.beta1),floor(num_outlier1 * options.beta2)); % estimation of the outlier number
    else  % same score, increase the estimation (num_outlier1) but can not go over num_outlier2
          % lower score, 
        num_outlier1    = min(floor(num_outlier1 * options.beta2),floor(num_outlier2 * 1)); % estimation of the outlier number        
    end
        
    num_outlier(iter+1,:)    = [num_outlier1, num_outlier2, num_outlier3, num_outlier4];

    if num_outlier2 == num_outlier1
        if num_outlier2 <= num_outlier_B
            diff_top    = index(num_outlier1 + 1:numsample); % 把top p%的outlier去掉
            data2       = data(diff_top,:);
            if options.adaHodge
                score       = AdaHodgerank(data2,score)'; %去掉top p%的outlier后再算L2 score
            else
                score       = Hodgerank(data2,score)'; %去掉top p%的outlier后再算L2 score
            end
            diff         = max(score(data(:,2))' - score(data(:,1))', 0);
            [diff_sort index]   = sort(diff,'descend');
            outlier_detect(:,iter+1) = diff > diff_sort(num_outlier1+1);        
        else
            if options.adaHodge
                [temp, index] = sort(score_B,'descend');
                score(index) = [n-1:-1:0];
            else
                score = score_B;
            end
            diff         = max(score(data(:,2))' - score(data(:,1))', 0);
            [diff_sort index]   = sort(diff,'descend');
            outlier_detect(:,iter+1) = diff > diff_sort(num_outlier1+1);        
        end
        break
    end
    
    if iter == options.max_iter
        break
    end

    score = score - min(score);
    score = (n - 1) * score / max(score);
    [score_sort, index2] = sort(score,'descend');
    score(index2(1)) = n - 1;
    delay_total = 0;
    delay = 0;
    for k = 2:n
        if score_sort(k-1) - score_sort(k) > options.intercept
            score(index2(k)) = score(index2(k-1)) - 1 - delay;
            delay = 0;
        else
            score(index2(k)) = score(index2(k-1));            
            delay = delay + 1;
            delay_total = delay_total + 1;
        end
    end
end

[score_sort, index2] = sort(score,'descend');

output.adjust = [];
for i = 2:length(score)
    if Z0(index2(i-1),index2(i)) < Z0(index2(i),index2(i-1))
        temp = score(index2(i));
        score(index2(i)) = score(index2(i-1));
        score(index2(i-1)) = temp;
        output.adjust = [output.adjust; index2(i-1), index2(i)];
    end
end

if ~isempty(output.adjust)
diff           = max(score(data(:,2))' - score(data(:,1))', 0);
data2 = data2(diff>0,:);
if options.adaHodge
    score       = AdaHodgerank(data2,score)'; %去掉top p%的outlier后再算L2 score
else
    score       = Hodgerank(data2,score)'; %去掉top p%的outlier后再算L2 score
end

end
outlier_detect = [outlier_detect diff > 0];   

output.Z0             = Z0;
output.num_outlier    = num_outlier;
output.outlier_detect = outlier_detect;
end