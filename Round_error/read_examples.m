%% READ_EXAMPLES :: read a given example files ('name'+'_i,j,s,c,p.dat') to produce 
%  The needed information to solve it (F,I,J,G,n,d,k)
%  The tag 'MIN' or 'MAX' slect the choice of F for the minimization or the
%  maximization problem
function [ F,I,J ,G ,n,d,k ] = read_examples(name,tag)

%% powers
s = strcat('Examples/',name,'/',name,'_p.dat');
p = dlmread(s);
%% d = max degree ;  k=0
d = max(sum(p'));
k=0;
%% coefficients
s = strcat('Examples/',name,'/',name,'_c.dat');
c = dlmread(s);
%% building F
if(strcmp(tag,'MAX'))
    F = [p -c']; % max(f(x)) =  -min(-f(x))
else
    F = [p c'];
end    


%% other infos
s = strcat('Examples/',name,'/',name,'_s.dat');
info = dlmread(s);
n = info(1);
complex_sparse = info(2);
n_semialg = info(3);

%% I
if(complex_sparse)
    %% Complex sparcisity pattern for I
    s = strcat('Examples/',name,'/',name,'_ic.dat');
    It = dlmread(s);
    I = cell(1,size(It,1));
    for j=1:size(It,1)
        index  = find(It(j,:)==0,1,'first');
        if(index>1)
            I{j} = It(j,1:(index-1));
        else
             I{j} = It(j,:);
        end
    end
    %% Complex sparsicity pattern for J (box edition !)
    J = I;
else
    s = strcat('Examples/',name,'/',name,'_i.dat');
    i = dlmread(s);
    for m=1 : size(i,1);
        I(m,:) = {i(m,:)};
    end
    s = strcat('Examples/',name,'/',name,'_j.dat');
    j = dlmread(s);
    for m=1 : size(j,1);
        J(m,:) = {j(m,:)};
    end
end

%% G
s = strcat('Examples/',name,'/',name,'_g.dat');
g = dlmread(s);
for m=1 : size(g,1);
    G(m,:) = {g(m,:)};
end

fin = size(g,1);

%% Rebuilding the semi-algebraic constraints
for m=1:n_semialg
   g = [];
   s = strcat('Examples/',name,'/',name,'_p_Gsa_',num2str(m),'.dat');
   ploc = dlmread(s);
   s = strcat('Examples/',name,'/',name,'_c_Gsa_',num2str(m),'.dat');
   cloc = dlmread(s) ; 
   for nm=1:length(cloc)
      lg = [ploc(nm,:) cloc(nm)]
      g = [g;lg]
   end  
   G(fin+m,:) =  {g};
end    

G


end

