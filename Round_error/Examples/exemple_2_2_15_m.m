%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 2;
name = 'exemple_2_2_15';
nparam = 19;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =(( 32 * x3 + ( 32 * x5 + ( 16 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11 + ( 7 * x12 + ( 8 * x13 + ( 9 * x14 + ( 10 * x15 + ( 11 * x16 + ( 12 * x17 + ( 13 * x18 + ( 14 * x19 + ( 15 * x20 + ( 16 * x21)))))))))))))))))) * x1 + (( 32 * x3 + ( 32 * x4 + ( 64 * x5 + ( 32 * x6 + ( 4 * x7 + ( 6 * x8 + ( 8 * x9 + ( 10 * x10 + ( 12 * x11 + ( 14 * x12 + ( 16 * x13 + ( 18 * x14 + ( 20 * x15 + ( 22 * x16 + ( 24 * x17 + ( 26 * x18 + ( 28 * x19 + ( 30 * x20 + ( 32 * x21))))))))))))))))))) * x2)) * x1 + (( 32 * x4 + ( 32 * x5 + ( 16 * x6 + ( 2 * x7 + ( 3 * x8 + ( 4 * x9 + ( 5 * x10 + ( 6 * x11 + ( 7 * x12 + ( 8 * x13 + ( 9 * x14 + ( 10 * x15 + ( 11 * x16 + ( 12 * x17 + ( 13 * x18 + ( 14 * x19 + ( 15 * x20 + ( 16 * x21)))))))))))))))))) * x2^2);
interval=[repmat([-1 1],nvar,1);repmat([-1 1],nparam,1)];
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


