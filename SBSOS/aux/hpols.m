%%***********************************************************%%
% This function generates all polynomials h_{\alpha \beta}    %
% described in Section The Sparse Bounded-SOS-hierarchy       %
%%***********************************************************%%

function [HABd,nAB] = hpols(G,d)
m = length(G);
np1 = size(G{1},2);
n = np1-1;
pows = cell(2*m,1);

for i = 1:m
    pows{i} = powpol(G{i},d,2);
    Gpow = G{i}(:,1:n);
    Gcof = G{i}(:,np1);
    Gpow = [zeros(1,n);Gpow];
    Gcof = [1;-Gcof];
    GG = [Gpow,Gcof];
    pows{i+m} = powpol(GG,d,2);
end

AB = monomials(2*m,d);
nAB = size(AB, 1);
ABp1 = AB + ones(nAB,2*m);
%% ANCIENNE VERSION
%  HABd = cell(nAB,1);
%%

%% AJOUT VECTORIZATION ALEXANDRE ROCCA
i = ABp1(:,1);
HABd = pows{1}(i(:,1));
H = HABd;
E = ABp1(1:nAB,2:2*m)>1; %% cropping on 2:2*m done here

[l,c] = find(E>0);
I = [l c];

for j = 2: 2*m
    J = find(I(:,2)==j-1); %% j-1 because of the above cropping
    if(size(J)>0)    
        HABd(I(J,1),1) = cellfun(@multpol, HABd(I(J,1),1), pows{j}(ABp1(I(J,1),j)), 'UniformOutput', 0);
    end
end
%%

%% ANCIENNE VERSION
% for i = 1: nAB
%     h = pows{1}{ABp1(i,1)};
%     for j = 2: 2*m
%         e = ABp1(i,j);
%         if e > 1
%            h = multpol(h,pows{j}{e});
%         end
%     end
%     HABd{i}=h;
% end
%%





