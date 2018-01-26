

%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 1;
name = 'sqroot';
nparam = 15;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q = (((( (-5/32) * x2 + ( (-1/32) * x5 + ( (-5/128) * x12 + ( (-5/128) * x13 + ( (-5/128) * x14 + ( (-5/128) * x15 + ( (-5/128) * x16))))))) * x1 + ( (3/16) * x2 + ( (1/16) * x4 + ( (1/16) * x8 + ( (1/16) * x9 + ( (1/16) * x10 + ( (1/16) * x11 + ( (1/16) * x16)))))))) * x1 + ( (-1/4) * x2 + ( (-1/8) * x3 + ( (-1/8) * x5 + ( (-1/8) * x6 + ( (-1/8) * x7 + ( (-1/8) * x11 + ( (-1/8) * x16)))))))) * x1 + (x2 + ( (1/2) * x3 + ( (3/2) * x4 + ( (3/2) * x7 + ( (3/2) * x11 + ( (3/2) * x16))))))) * x1 + (x4 + (x7 + (x11 + (x16))));

interval=[repmat([0 1],nvar,1);repmat([-1 1],nparam,1)];
qsdp = box_norm(q,vars,interval);
%% Computation of the power matrix and the coefficient list
[powers,coefficients] = getexponentbase(qsdp,vars);
p = str2double(sdisplay(powers));
c = str2double(sdisplay(coefficients));

n = nvar+nparam;
[I,J] = build_box_sparcity(nvar,nparam);
cstr = []
n_semialg =  size(cstr,1);
system_info = [n complex_sparse n_semialg];
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


