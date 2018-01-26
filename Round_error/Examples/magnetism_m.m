

%% Buildind the model files for schwefel model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 7;
name = 'magnetism';
nparam = 27;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =(( (2/1) * x8 + (x15 + (x18 + (x21 + (x24 + (x27 + (x30 + (x33 + (x34))))))))) * x1 + ( - x8 + ( - x34))) * x1 + (( (4/1) * x9 + ( (2/1) * x16 + ( (2/1) * x17 + ( (2/1) * x18 + ( (2/1) * x21 + ( (2/1) * x24 + ( (2/1) * x27 + ( (2/1) * x30 + ( (2/1) * x33 + ( (2/1) * x34)))))))))) * x2^2 + (( (4/1) * x10 + ( (2/1) * x19 + ( (2/1) * x20 + ( (2/1) * x21 + ( (2/1) * x24 + ( (2/1) * x27 + ( (2/1) * x30 + ( (2/1) * x33 + ( (2/1) * x34))))))))) * x3^2 + (( (4/1) * x11 + ( (2/1) * x22 + ( (2/1) * x23 + ( (2/1) * x24 + ( (2/1) * x27 + ( (2/1) * x30 + ( (2/1) * x33 + ( (2/1) * x34)))))))) * x4^2 + (( (4/1) * x12 + ( (2/1) * x25 + ( (2/1) * x26 + ( (2/1) * x27 + ( (2/1) * x30 + ( (2/1) * x33 + ( (2/1) * x34))))))) * x5^2 + (( (4/1) * x13 + ( (2/1) * x28 + ( (2/1) * x29 + ( (2/1) * x30 + ( (2/1) * x33 + ( (2/1) * x34)))))) * x6^2 + (( (4/1) * x14 + ( (2/1) * x31 + ( (2/1) * x32 + ( (2/1) * x33 + ( (2/1) * x34))))) * x7^2))))));
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


