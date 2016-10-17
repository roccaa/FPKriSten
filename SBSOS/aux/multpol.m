%%***************************************************************%%
% This function builds the product r of to polynomials p and q.   %
% The resulting polynomial r is not simplified, i.e. r is likely  %
% to have several entries for the same monomial.                  %
%%***************************************************************%%

function r = multpol(p,q)

[mp, np1] = size(p);
[mq, nq1] = size(q);
if ~(np1 ==nq1)
    error('number of variables does not coincide')
end

%% Modification Alexandre ROCCA
% r = zeros(mp*mq,np1);
%%
n = np1-1;
P = p(:,1:n);
Q = q(:,1:n);

%% AJOUT VECTORIZATION ALEXANDRE ROCCA
r = repmat(Q,mp,1);
r3 = repmat(q(:,np1),mp,1);
Prep = P(repmat(1:size(p,1),mq,1),:) ; 
r = r+Prep;
t = p(:,np1);
prep = t(repmat(1:size(p,1),mq,1),:) ; 
r3 = r3.*prep;
r = [r r3];
%% ANCIENE VERSION
% for i = 1:mp
%     for j = 1:mq
%         r((i-1)*mq +j,1:n) = P(i,:)+Q(j,:);
%         v((i-1)*mq +j,1:n) = Q(j,:);
%         w((i-1)*mq +j,1:n) = P(i,:);
%         r((i-1)*mq +j,np1) = p(i,np1)*q(j,np1);
%     end
% end



end
