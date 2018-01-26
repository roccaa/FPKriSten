

%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 3;
name = 'rigidBody2';
nparam = 15;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q = (((( - x4 + ( (-2/1) * x5 + ( - x6 + ( - x13 + ( - x14 + ( - x15 + ( - x16 + ( - x17 + ( - x18))))))))) * x3) * x2 + (( (2/1) * x4 + ( (2/1) * x5 + ( (2/1) * x6 + ( (2/1) * x7 + ( (2/1) * x8 + ( (2/1) * x9 + ( (2/1) * x12 + ( (2/1) * x16 + ( (2/1) * x17 + ( (2/1) * x18)))))))))) * x3)) * x2) * x1 + (( - x5 + ( - x18)) * x2 + (( (12/1) * x6 + ( (6/1) * x10 + ( (6/1) * x11 + ( (3/1) * x12 + ( (3/1) * x16 + ( (6/1) * x17 + ( (6/1) * x18))))))) * x3^2));
interval=[repmat([-15 15],nvar,1);repmat([-1 1],nparam,1)];
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


