%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 2;
name = 'exemple_2_2_20';
nparam = 24;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =(( (42/1) * x3 + ( (42/1) * x5 + ( (21/1) * x6 + ( (2/1) * x7 + ( (3/1) * x8 + ( (4/1) * x9 + ( (5/1) * x10 + ( (6/1) * x11 + ( (7/1) * x12 + ( (8/1) * x13 + ( (9/1) * x14 + ( (10/1) * x15 + ( (11/1) * x16 + ( (12/1) * x17 + ( (13/1) * x18 + ( (14/1) * x19 + ( (15/1) * x20 + ( (16/1) * x21 + ( (17/1) * x22 + ( (18/1) * x23 + ( (19/1) * x24 + ( (20/1) * x25 + ( (21/1) * x26))))))))))))))))))))))) * x1 + (( (42/1) * x3 + ( (42/1) * x4 + ( (84/1) * x5 + ( (42/1) * x6 + ( (4/1) * x7 + ( (6/1) * x8 + ( (8/1) * x9 + ( (10/1) * x10 + ( (12/1) * x11 + ( (14/1) * x12 + ( (16/1) * x13 + ( (18/1) * x14 + ( (20/1) * x15 + ( (22/1) * x16 + ( (24/1) * x17 + ( (26/1) * x18 + ( (28/1) * x19 + ( (30/1) * x20 + ( (32/1) * x21 + ( (34/1) * x22 + ( (36/1) * x23 + ( (38/1) * x24 + ( (40/1) * x25 + ( (42/1) * x26)))))))))))))))))))))))) * x2)) * x1 + (( (42/1) * x4 + ( (42/1) * x5 + ( (21/1) * x6 + ( (2/1) * x7 + ( (3/1) * x8 + ( (4/1) * x9 + ( (5/1) * x10 + ( (6/1) * x11 + ( (7/1) * x12 + ( (8/1) * x13 + ( (9/1) * x14 + ( (10/1) * x15 + ( (11/1) * x16 + ( (12/1) * x17 + ( (13/1) * x18 + ( (14/1) * x19 + ( (15/1) * x20 + ( (16/1) * x21 + ( (17/1) * x22 + ( (18/1) * x23 + ( (19/1) * x24 + ( (20/1) * x25 + ( (21/1) * x26))))))))))))))))))))))) * x2^2);
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


