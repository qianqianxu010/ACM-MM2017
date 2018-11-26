function score = AdaHodgerank(incomp, score, alpha, delta)
%--------------------------------------------------------------------------
% AdaHodgerank.m: Adaptive HodgeRank by iteratively updating the score
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Adaptive HodgeRank by iteratively updating the score
%
% USAGE: 
%    score = AdaHodgerank(incomp, score, alpha, delta)
%    score = AdaHodgerank(incomp, score)
%    score = AdaHodgerank(incomp)
%
% INPUT ARGUMENTS:
% incomp      Input matrix with paired comparison, the matrix has 
%             2 columns, and for each row, the rank of the first column
%             is higher than the second column for this comparison.
% score       Initial Score, default is all ones.
% alpha       The initial gap between the score to be consider
%             has different score.
% delta       The change of alpha for each iteration.
%
% OUTPUT ARGUMENTS:
% score       The final score which are integers from 0 to n-1 
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
% Hodgerank, AODHodgerank
%
%%
n       = max(incomp(:));

if nargin < 2
    score = ones(n,1);
end
if nargin < 3
    alpha = 0.5;
end
if nargin < 4
    delta = 0.9;
end

for iter = 1:300
    if max(score) == min(score)
        delay_total = n-1;
        score = ones(n,1);
    else
        score = score - min(score);
        score = (n-1) * score / max(score);
        [score_sort, index] = sort(score,'descend');
        score(index(1)) = n-1;
        delay_total = 0;
        delay = 0;
        for i = 2:n
            if score_sort(i-1) - score_sort(i) > alpha
                score(index(i)) = score(index(i-1)) - 1 - delay;
                delay = 0;
            else
                score(index(i)) = score(index(i-1));            
                delay = delay + 1;
                delay_total = delay_total + 1;
            end
        end
    end
    score = Hodgerank(incomp, score);    
    diff  = max(score(incomp(:,2))' - score(incomp(:,1))', 0);
    if sum(diff) == 0
       [temp, index] = sort(score,'descend');
       score(index) = [n-1:-1:0];
       break;
    end

    alpha = delta * alpha;
    if delay_total == 0
        break
    end
end
score = score - min(score);