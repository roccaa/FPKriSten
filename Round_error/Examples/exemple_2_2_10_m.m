%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 2;
name = 'exemple_2_2_10';
nparam = 14;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =(( 22 * x3 + ( 22 * x5 + ( 11 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11 + ( 7 * x12 + ( 8 * x13 + ( 9 * x14 + ( 10 * x15 + ( 11 * x16))))))))))))) * x1 + (( 22 * x3 + ( 22 * x4 + ( 44 * x5 + ( 22 * x6 + ( 4 * x7 + ( 6 * x8 + ( 8 * x9 + ( 10 * x10 + ( 12 * x11 + ( 14 * x12 + ( 16 * x13 + ( 18 * x14 + ( 20 * x15 + ( 22 * x16)))))))))))))) * x2)) * x1 + (( 22 * x4 + ( 22 * x5 + ( 11 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11 + ( 7 * x12 + ( 8 * x13 + ( 9 * x14 + ( 10 * x15 + ( 11 * x16))))))))))))) * x2^2);
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


