%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 10;
name = 'exemple_10_2_2';
nparam = 22;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q =(( 6 * x11 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32))))))))))))) * x1 + (( 6 * x11 + ( 6 * x12 + ( 12 * x21 + ( 12 * x22 + ( 12 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x2 + (( 6 * x11 + ( 6 * x13 + ( 6 * x21 + ( 12 * x22 + ( 12 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x3 + (( 6 * x11 + ( 6 * x14 + ( 6 * x21 + ( 6 * x22 + ( 12 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x4 + (( 6 * x11 + ( 6 * x15 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x5 + (( 6 * x11 + ( 6 * x16 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x6 + (( 6 * x11 + ( 6 * x17 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x7 + (( 6 * x11 + ( 6 * x18 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x8 + (( 6 * x11 + ( 6 * x19 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x9 + (( 6 * x11 + ( 6 * x20 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x10)))))))))) * x1 + ((( 6 * x12 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32))))))))))))) * x2 + (( 6 * x12 + ( 6 * x13 + ( 6 * x21 + ( 12 * x22 + ( 12 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x3 + (( 6 * x12 + ( 6 * x14 + ( 6 * x21 + ( 6 * x22 + ( 12 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x4 + (( 6 * x12 + ( 6 * x15 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x5 + (( 6 * x12 + ( 6 * x16 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x6 + (( 6 * x12 + ( 6 * x17 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x7 + (( 6 * x12 + ( 6 * x18 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x8 + (( 6 * x12 + ( 6 * x19 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x9 + (( 6 * x12 + ( 6 * x20 + ( 6 * x21 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))))) * x10))))))))) * x2 + ((( 6 * x13 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32)))))))))))) * x3 + (( 6 * x13 + ( 6 * x14 + ( 6 * x22 + ( 12 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))))) * x4 + (( 6 * x13 + ( 6 * x15 + ( 6 * x22 + ( 6 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))))) * x5 + (( 6 * x13 + ( 6 * x16 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))))) * x6 + (( 6 * x13 + ( 6 * x17 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))))) * x7 + (( 6 * x13 + ( 6 * x18 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))))) * x8 + (( 6 * x13 + ( 6 * x19 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))))) * x9 + (( 6 * x13 + ( 6 * x20 + ( 6 * x22 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))))) * x10)))))))) * x3 + ((( 6 * x14 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32))))))))))) * x4 + (( 6 * x14 + ( 6 * x15 + ( 6 * x23 + ( 12 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))) * x5 + (( 6 * x14 + ( 6 * x16 + ( 6 * x23 + ( 6 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))) * x6 + (( 6 * x14 + ( 6 * x17 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))) * x7 + (( 6 * x14 + ( 6 * x18 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))) * x8 + (( 6 * x14 + ( 6 * x19 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))) * x9 + (( 6 * x14 + ( 6 * x20 + ( 6 * x23 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))))) * x10))))))) * x4 + ((( 6 * x15 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32)))))))))) * x5 + (( 6 * x15 + ( 6 * x16 + ( 6 * x24 + ( 12 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))) * x6 + (( 6 * x15 + ( 6 * x17 + ( 6 * x24 + ( 6 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))) * x7 + (( 6 * x15 + ( 6 * x18 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))) * x8 + (( 6 * x15 + ( 6 * x19 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))) * x9 + (( 6 * x15 + ( 6 * x20 + ( 6 * x24 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))))) * x10)))))) * x5 + ((( 6 * x16 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32))))))))) * x6 + (( 6 * x16 + ( 6 * x17 + ( 6 * x25 + ( 12 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))) * x7 + (( 6 * x16 + ( 6 * x18 + ( 6 * x25 + ( 6 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))) * x8 + (( 6 * x16 + ( 6 * x19 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))) * x9 + (( 6 * x16 + ( 6 * x20 + ( 6 * x25 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))))) * x10))))) * x6 + ((( 6 * x17 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32)))))))) * x7 + (( 6 * x17 + ( 6 * x18 + ( 6 * x26 + ( 12 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))) * x8 + (( 6 * x17 + ( 6 * x19 + ( 6 * x26 + ( 6 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))) * x9 + (( 6 * x17 + ( 6 * x20 + ( 6 * x26 + ( 6 * x27 + ( 6 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))))) * x10)))) * x7 + ((( 6 * x18 + ( 6 * x27 + ( 6 * x28 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32))))))) * x8 + (( 6 * x18 + ( 6 * x19 + ( 6 * x27 + ( 12 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))) * x9 + (( 6 * x18 + ( 6 * x20 + ( 6 * x27 + ( 6 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32)))))))) * x10))) * x8 + ((( 6 * x19 + ( 6 * x28 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32)))))) * x9 + (( 6 * x19 + ( 6 * x20 + ( 6 * x28 + ( 12 * x29 + ( 6 * x30 + ( 4 * x31 + ( 6 * x32))))))) * x10)) * x9 + (( 6 * x20 + ( 6 * x29 + ( 3 * x30 + ( 2 * x31 + ( 3 * x32))))) * x10^2)))))))));
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


