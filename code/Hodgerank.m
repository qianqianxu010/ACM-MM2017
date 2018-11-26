function [score, Z] = Hodgerank(data, score)

%--------------------------------------------------------------------------
% Hodgerank.m: HodgeRank
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    HodgeRank
%
% USAGE: 
%    score = Hodgerank(data, score)
%    score = Hodgerank(data)
%
% INPUT ARGUMENTS:
% data        Input matrix with paired comparison, the matrix has 
%             2 columns, and for each row, the rank of the first column
%             is higher than the second column for this comparison.
% score       Initial Score, default is all ones.
%
% OUTPUT ARGUMENTS:
% score       The output score
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
% AdaHodgerank, AODHodgerank
%
%%
n       = max(data(:)); % ID is 1 to n

if nargin < 2
    score = ones(n,1);
end

[compare,temp] = size(data); % compare is the number of comparison

% transform to (d,y,w)
Z = zeros(n,n);
for k = 1:compare
    a = data(k,:);
    Z(a(1),a(2)) = Z(a(1),a(2)) + 1;
end

m = sum(sum(Z~=0));
d = zeros(m,n);
w = zeros(m,1);
y = ones(m,1);

k = 0;
for i = 1:n
    for j = 1:n
        if (i~=j && Z(i,j)~=0)
            k = k+1;
            w(k) = Z(i,j);
            d(k,i) = 1;
            d(k,j) = -1;
            y(k)  =  max(score(i)-score(j),1);
        end
    end
end
display(['m is ' num2str(m), ', n(n-1)/2 should be ' num2str(n*(n-1)/2)])
score  = lsqr(d'*diag(w)*d,d'*(w.*y));
%score  = score - min(score);