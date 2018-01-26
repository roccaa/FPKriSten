%% Buildind the model files for floudas4_7 model
%% Normalization (from box representation) in [0:1]
%% Use Yalmip
use_semialgebraic = 1;
complex_sparse = 0;
nvar = 2;
name ='floudas4_7';
nparam = 8;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);

% Model
q = ( -12 * x3 + ( -12 * x5 + ( -12 * x6 + ( -12 * x8 + ( -12 * x10))))) * x1 + ((( 2 * x4 + (x9 + (x10))) * x2 + ( -7 * x4 + ( -7 * x7 + ( -7 * x8 + ( -7 * x10))))) * x2);
cstr = [(-2*x1^4+2-x2)/2];

n_semialg =  length(cstr);
    
interval=[[0 2];[0 3];repmat([-1 1],nparam,1)];
qsdp = box_norm(q,vars,interval);

% normalizing semi-algebraic equation
for i=1:n_semialg
    cstr(i) = box_norm(cstr(i),vars,interval);
end
%% Computation of the power matrix and the coefficient list
% qsdp = eval(q_norm);
[powers,coefficients] = getexponentbase(qsdp,vars);
p = str2double(sdisplay(powers));
c = str2double(sdisplay(coefficients));



n = nvar+nparam;
[I,J] = build_box_sparcity(nvar,nparam);
% Completing J with semi-algebraic constraints
ci = nvar+nparam+1:nvar+nparam+n_semialg;

C = repmat(ci,nparam,1);
J = [J C];

system_info = [n complex_sparse n_semialg];
G = create_unitBox(n);

%% Similar to box only
path = [name '/' name];
mkdir(name);
dlmwrite([path '_s.dat'],system_info);
dlmwrite([path '_c.dat'],c);
dlmwrite([path '_p.dat'],p);
%% Classical pattern
dlmwrite([path '_i.dat'],I);
dlmwrite([path '_j.dat'],J);


%% getting info on Semialgebraic
for i=1:n_semialg
    [powersloc,coefficientsloc] = getexponentbase(cstr(i),vars);
     ploc = str2double(sdisplay(powersloc));
     cloc = str2double(sdisplay(coefficientsloc));
     dlmwrite([path '_p_Gsa_' num2str(i) '.dat'],ploc);
     dlmwrite([path '_c_Gsa_' num2str(i) '.dat'],cloc);
end

dlmwrite([path '_g.dat'],G);


