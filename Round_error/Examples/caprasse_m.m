

%% Buildind the model files for schwefel model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 4;
name = 'caprasse';
nparam = 34;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =((( - x5 + ( (-3/1) * x7 + ( - x9 + ( - x10 + ( - x11 + ( - x12 + ( - x17 + ( - x22 + ( - x27 + ( - x28 + ( - x31 + ( - x34 + ( - x37 + ( - x38)))))))))))))) * x3^2 + (( (4/1) * x5 + ( (4/1) * x7 + ( (8/1) * x8 + ( (4/1) * x18 + ( (4/1) * x19 + ( (4/1) * x20 + ( (4/1) * x21 + ( (4/1) * x22 + ( (4/1) * x27 + ( (4/1) * x28 + ( (4/1) * x31 + ( (4/1) * x34 + ( (4/1) * x37 + ( (4/1) * x38)))))))))))))) * x4^2 + ( (4/1) * x5 + ( (4/1) * x7 + ( (4/1) * x18 + ( (4/1) * x19 + ( (4/1) * x28 + ( (4/1) * x31 + ( (4/1) * x34 + ( (4/1) * x37 + ( (4/1) * x38))))))))))) * x3) * x1 + (((( (4/1) * x6 + ( (8/1) * x7 + ( (4/1) * x8 + ( (4/1) * x13 + ( (4/1) * x14 + ( (4/1) * x15 + ( (4/1) * x16 + ( (4/1) * x17 + ( (4/1) * x22 + ( (4/1) * x27 + ( (4/1) * x28 + ( (4/1) * x31 + ( (4/1) * x34 + ( (4/1) * x37 + ( (4/1) * x38))))))))))))))) * x4) * x3^2 + ((( (2/1) * x6 + ( (6/1) * x8 + ( (2/1) * x23 + ( (2/1) * x24 + ( (2/1) * x25 + ( (2/1) * x26 + ( (2/1) * x27 + ( (2/1) * x28 + ( (2/1) * x31 + ( (2/1) * x34 + ( (2/1) * x37 + ( (2/1) * x38)))))))))))) * x4^2 + ( (-10/1) * x6 + ( (-10/1) * x8 + ( (-10/1) * x32 + ( (-10/1) * x33 + ( (-10/1) * x34 + ( (-10/1) * x37 + ( (-10/1) * x38)))))))) * x4)) * x2 + (( (8/1) * x7 + ( (4/1) * x29 + ( (4/1) * x30 + ( (4/1) * x31 + ( (4/1) * x34 + ( (4/1) * x37 + ( (4/1) * x38))))))) * x3^2 + (( (-20/1) * x8 + ( (-10/1) * x35 + ( (-10/1) * x36 + ( (-10/1) * x37 + ( (-10/1) * x38))))) * x4^2 + ( (2/1) * x38))));
interval=[repmat([-0.5 0.5],nvar,1);repmat([-1 1],nparam,1)];
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


