%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 2;
name = 'exemple_2_2_5';
nparam = 9;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =(( 12 * x3 + ( 12 * x5 + ( 6 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11)))))))) * x1 + (( 12 * x3 + ( 12 * x4 + ( 24 * x5 + ( 12 * x6 + ( 4 * x7 + ( 6 * x8 + ( 8 * x9 + ( 10 * x10 + ( 12 * x11))))))))) * x2)) * x1 + (( 12 * x4 + ( 12 * x5 + ( 6 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11)))))))) * x2^2);
interval=[repmat([-1 1],nvar,1);repmat([-1 1],nparam,1)];
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


