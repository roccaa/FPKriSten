

%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 1;
name = 'sineOrder3';
nparam = 6;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q = ((( (-5120184131/10000000000) * x2 + ( (-1290061377/10000000000) * x4 + ( (-1290061377/10000000000) * x5 + ( (-1290061377/10000000000) * x6 + ( (-1290061377/10000000000) * x7))))) * x1 + ( (1/2) * x2 + ( (1/2) * x3 + ( (1/2) * x7)))) * x1 + ( (4774648293/5000000000) * x2 + ( (4774648293/5000000000) * x3 + ( (4774648293/5000000000) * x7)))) * x1;

interval=[repmat([-2 2],nvar,1);repmat([-1 1],nparam,1)];
qsdp = box_norm(q,vars,interval);
%% Computation of the power matrix and the coefficient list
[powers,coefficients] = getexponentbase(qsdp,vars);
p = str2double(sdisplay(powers));
c = str2double(sdisplay(coefficients));

n = nvar+nparam;
[I,J] = build_box_sparcity(nvar,nparam);
system_info = [n complex_sparse];
G = create_unitBox(n);

path = [name '/' name];
mkdir(name);
dlmwrite([path '_g.dat'],G);
dlmwrite([path '_s.dat'],system_info);
dlmwrite([path '_c.dat'],c);
dlmwrite([path '_p.dat'],p);

%% Classical pattern
dlmwrite([path '_i.dat'],I);
dlmwrite([path '_j.dat'],J);


