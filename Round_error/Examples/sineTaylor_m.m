

%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 1;
name = 'sineTaylor';
nparam = 13;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q = (((( (-1/720) * x2 + ( (-1/5040) * x3 + ( (-1/5040) * x4 + ( (-1/5040) * x7 + ( (-1/5040) * x8 + ( (-1/5040) * x11 + ( (-1/5040) * x12 + ( (-1/5040) * x13 + ( (-1/5040) * x14))))))))) * x1^2 + ( (1/24) * x2 + ( (1/120) * x3 + ( (1/120) * x4 + ( (1/120) * x7 + ( (1/120) * x8 + ( (1/120) * x9 + ( (1/120) * x10 + ( (1/120) * x14))))))))) * x1^2 + ( (-1/2) * x2 + ( (-1/6) * x3 + ( (-1/6) * x4 + ( (-1/6) * x5 + ( (-1/6) * x6 + ( (-1/6) * x10 + ( (-1/6) * x14)))))))) * x1^2 + (x2 + (x6 + (x10 + (x14))))) * x1

interval=[repmat([-1.57079632679 1.57079632679],nvar,1);repmat([-1 1],nparam,1)];
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


