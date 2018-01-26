
%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 6;
name = 'kepler2';
nparam = 42;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q = ((( (-2/1) * x7 + ( - x10 + ( - x13 + ( - x14 + ( - x15 + ( - x16 + ( - x17 + ( - x18 + ( - x19 + ( - x20 + ( - x28 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48)))))))))))))))) * x4) * x1 + (((x7 + (x8 + (x10 + (x13 + (x15 + (x16 + (x17 + (x18 + (x19 + (x20 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))))))) * x4 + ((x7 + (x8 + (x11 + (x21 + (x22 + (x23 + (x24 + (x25 + (x26 + (x27 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))))))) * x5 + (( - x7 + ( - x8 + ( - x12 + ( - x43 + ( - x44 + ( - x45 + ( - x48))))))) * x6))) * x2 + (((x7 + (x9 + (x10 + (x13 + (x16 + (x17 + (x18 + (x19 + (x20 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x4 + (( - x7 + ( - x9 + ( - x11 + ( - x40 + ( - x41 + ( - x42 + ( - x45 + ( - x48)))))))) * x5 + ((x7 + (x9 + (x12 + (x29 + (x30 + (x31 + (x32 + (x33 + (x34 + (x35 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x6))) * x3 + ((( - x7 + ( (-2/1) * x10 + ( - x13 + ( - x17 + ( - x18 + ( - x19 + ( - x20 + ( - x28 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48))))))))))))) * x4 + ((x7 + (x10 + (x11 + (x13 + (x18 + (x19 + (x20 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))) * x5 + ((x7 + (x10 + (x12 + (x13 + (x19 + (x20 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))) * x6))) * x4)))) * x1 + (((( (-2/1) * x8 + ( - x11 + ( - x21 + ( - x22 + ( - x23 + ( - x24 + ( - x25 + ( - x26 + ( - x27 + ( - x28 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48))))))))))))))) * x5) * x2 + ((( - x8 + ( - x9 + ( - x10 + ( - x37 + ( - x38 + ( - x39 + ( - x42 + ( - x45 + ( - x48))))))))) * x4 + ((x8 + (x9 + (x11 + (x21 + (x23 + (x24 + (x25 + (x26 + (x27 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x5 + ((x8 + (x9 + (x12 + (x29 + (x30 + (x31 + (x32 + (x33 + (x34 + (x35 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x6))) * x3 + (((x8 + (x10 + (x11 + (x21 + (x24 + (x25 + (x26 + (x27 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))))) * x5) * x4 + ((( - x8 + ( (-2/1) * x11 + ( - x21 + ( - x25 + ( - x26 + ( - x27 + ( - x28 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48)))))))))))) * x5 + ((x8 + (x11 + (x12 + (x21 + (x26 + (x27 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))) * x6)) * x5)))) * x2 + (((( (-2/1) * x9 + ( - x12 + ( - x29 + ( - x31 + ( - x32 + ( - x33 + ( - x34 + ( - x35 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48))))))))))))) * x6) * x3 + (((x9 + (x10 + (x12 + (x29 + (x32 + (x33 + (x34 + (x35 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))) * x6) * x4 + (((x9 + (x11 + (x12 + (x29 + (x33 + (x34 + (x35 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))) * x6) * x5 + (( - x9 + ( (-2/1) * x12 + ( - x29 + ( - x34 + ( - x35 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48)))))))))) * x6^2)))) * x3 + (((( - x10 + ( - x11 + ( - x12 + ( - x46 + ( - x47 + ( - x48)))))) * x6) * x5) * x4)));
interval=[repmat([4 6.36],nvar,1);repmat([-1 1],nparam,1)];
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


% %% I (&J) complex sparcisity pattern
% pattern = {[1,2,4];[1,3,4];[2,3,4]};
% P = nvar:(nvar+nparam);
% Crep=[];
% Crep = cell(size(pattern,1)*size(P,2),1);
% for i=1:size(P,2)
%     for j=1:size(pattern,1)
%         Crep{(i-1)*size(pattern,1) +j} = [pattern{j} P(i)];
%     end
% end


%% I (complex sparsisity pattern)
% dlmwrite([path '_ic.dat'],'');
% for i=1:size(Crep)
%     dlmwrite([path '_ic.dat']',Crep{i},'-append');
% end
%% J = I

%% Classical pattern
dlmwrite([path '_i.dat'],I);
dlmwrite([path '_j.dat'],J);


