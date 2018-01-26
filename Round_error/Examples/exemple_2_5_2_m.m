%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 2;
name = 'exemple_2_5_2';
nparam = 9;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =((((( 15 * x3 + ( 15 * x5 + ( 3 * x6 + ( 3 * x7 + ( 3 * x8 + ( 3 * x9 + ( 2 * x10 + ( 3 * x11)))))))) * x1 + (( 60 * x3 + ( 15 * x4 + ( 75 * x5 + ( 15 * x6 + ( 15 * x7 + ( 15 * x8 + ( 15 * x9 + ( 10 * x10 + ( 15 * x11))))))))) * x2)) * x1 + (( 90 * x3 + ( 60 * x4 + ( 150 * x5 + ( 30 * x6 + ( 30 * x7 + ( 30 * x8 + ( 30 * x9 + ( 20 * x10 + ( 30 * x11))))))))) * x2^2)) * x1 + (( 60 * x3 + ( 90 * x4 + ( 150 * x5 + ( 30 * x6 + ( 30 * x7 + ( 30 * x8 + ( 30 * x9 + ( 20 * x10 + ( 30 * x11))))))))) * x2^3)) * x1 + (( 15 * x3 + ( 60 * x4 + ( 75 * x5 + ( 15 * x6 + ( 15 * x7 + ( 15 * x8 + ( 15 * x9 + ( 10 * x10 + ( 15 * x11))))))))) * x2^4)) * x1 + (( 15 * x4 + ( 15 * x5 + ( 3 * x6 + ( 3 * x7 + ( 3 * x8 + ( 3 * x9 + ( 2 * x10 + ( 3 * x11)))))))) * x2^5);
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


