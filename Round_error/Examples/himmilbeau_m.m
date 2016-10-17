

%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 2;
name = 'himmilbeau';
nparam = 11;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =((( (4/1) * x3 + ( (2/1) * x5 + ( (2/1) * x6 + ( (2/1) * x7 + (x8 + (x13)))))) * x1^2 + (( (4/1) * x3 + ( (2/1) * x4 + ( (2/1) * x5 + ( (4/1) * x6 + ( (4/1) * x7 + ( (2/1) * x8 + ( (2/1) * x13))))))) * x2 + ( (-42/1) * x3 + ( (-22/1) * x5 + ( (-22/1) * x6 + ( (-44/1) * x7 + ( (-22/1) * x8 + ( (2/1) * x10 + ( (2/1) * x11 + (x12 + ( (-21/1) * x13))))))))))) * x1 + (( (2/1) * x3 + ( (4/1) * x4 + ( (2/1) * x9 + ( (4/1) * x10 + ( (4/1) * x11 + ( (2/1) * x12 + ( (2/1) * x13))))))) * x2^2 + ( (-14/1) * x3 + ( (-14/1) * x10 + ( (-28/1) * x11 + ( (-14/1) * x12 + ( (-14/1) * x13))))))) * x1 + (((( (4/1) * x4 + ( (2/1) * x9 + ( (2/1) * x10 + ( (2/1) * x11 + (x12 + (x13)))))) * x2^2 + ( (-26/1) * x4 + ( (2/1) * x6 + ( (2/1) * x7 + (x8 + ( (-14/1) * x9 + ( (-14/1) * x10 + ( (-28/1) * x11 + ( (-14/1) * x12 + ( (-13/1) * x13)))))))))) * x2 + ( (-22/1) * x4 + ( (-22/1) * x6 + ( (-44/1) * x7 + ( (-22/1) * x8 + ( (-22/1) * x13)))))) * x2 + ( (242/1) * x7 + ( (121/1) * x8 + ( (98/1) * x11 + ( (49/1) * x12 + ( (170/1) * x13))))));
interval=[repmat([-5 5],nvar,1);repmat([-1 1],nparam,1)];
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


