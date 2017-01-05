%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 5;
name = 'exemple_5_2_2';
nparam = 12;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =(( 6 * x6 + ( 6 * x11 + ( 6 * x12 + ( 6 * x13 + ( 6 * x14 + ( 3 * x15 + ( 2 * x16 + ( 3 * x17)))))))) * x1 + (( 6 * x6 + ( 6 * x7 + ( 12 * x11 + ( 12 * x12 + ( 12 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17))))))))) * x2 + (( 6 * x6 + ( 6 * x8 + ( 6 * x11 + ( 12 * x12 + ( 12 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17))))))))) * x3 + (( 6 * x6 + ( 6 * x9 + ( 6 * x11 + ( 6 * x12 + ( 12 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17))))))))) * x4 + (( 6 * x6 + ( 6 * x10 + ( 6 * x11 + ( 6 * x12 + ( 6 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17))))))))) * x5))))) * x1 + ((( 6 * x7 + ( 6 * x11 + ( 6 * x12 + ( 6 * x13 + ( 6 * x14 + ( 3 * x15 + ( 2 * x16 + ( 3 * x17)))))))) * x2 + (( 6 * x7 + ( 6 * x8 + ( 6 * x11 + ( 12 * x12 + ( 12 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17))))))))) * x3 + (( 6 * x7 + ( 6 * x9 + ( 6 * x11 + ( 6 * x12 + ( 12 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17))))))))) * x4 + (( 6 * x7 + ( 6 * x10 + ( 6 * x11 + ( 6 * x12 + ( 6 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17))))))))) * x5)))) * x2 + ((( 6 * x8 + ( 6 * x12 + ( 6 * x13 + ( 6 * x14 + ( 3 * x15 + ( 2 * x16 + ( 3 * x17))))))) * x3 + (( 6 * x8 + ( 6 * x9 + ( 6 * x12 + ( 12 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17)))))))) * x4 + (( 6 * x8 + ( 6 * x10 + ( 6 * x12 + ( 6 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17)))))))) * x5))) * x3 + ((( 6 * x9 + ( 6 * x13 + ( 6 * x14 + ( 3 * x15 + ( 2 * x16 + ( 3 * x17)))))) * x4 + (( 6 * x9 + ( 6 * x10 + ( 6 * x13 + ( 12 * x14 + ( 6 * x15 + ( 4 * x16 + ( 6 * x17))))))) * x5)) * x4 + (( 6 * x10 + ( 6 * x14 + ( 3 * x15 + ( 2 * x16 + ( 3 * x17))))) * x5^2))));
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


