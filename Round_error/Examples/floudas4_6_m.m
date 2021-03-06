  
    
%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
use_semialgebraic = 1;
complex_sparse = 0;
nvar = 2;
name ='floudas4_6';
nparam = 4;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);

% Model
q = ( - x3 + ( - x5 + ( - x6))) * x1 + (( - x4 + ( - x6)) * x2)
cstr = [(2*x1^4-8*x1^3+8*x1*x1-x2)/18;
		(4*x1^4-32*x1^3+88*x1*x1-96*x1+36-x2)/36]
    

n_semialg =  size(cstr,1);
    
interval=[[0 3];[0 4];repmat([-1 1],nparam,1)];
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


