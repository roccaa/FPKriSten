

%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 4;
name = 'kepler1';
nparam = 28;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);

q = ((( (-2/1) * x5 + ( - x8 + ( - x9 + ( - x10 + ( - x11 + ( - x12 + ( - x13 + ( - x14 + ( - x19 + ( - x24 + ( - x27 + ( - x29 + ( - x31 + ( - x32)))))))))))))) * x4) * x1 + (((x5 + (x6 + (x8 + (x9 + (x11 + (x12 + (x13 + (x14 + (x19 + (x24 + (x27 + (x29 + (x31 + (x32)))))))))))))) * x4 + (x15 + (x16 + (x17 + (x18 + (x19 + (x24 + (x27 + (x29 + ( - x30)))))))))) * x2 + (((x5 + (x7 + (x8 + (x9 + (x12 + (x13 + (x14 + (x19 + (x24 + (x27 + (x29 + (x31 + (x32))))))))))))) * x4 + (x20 + (x21 + (x22 + (x23 + (x24 + (x27 + ( - x28)))))))) * x3 + (( - x5 + ( (-2/1) * x8 + ( - x9 + ( - x13 + ( - x14 + ( - x19 + ( - x24 + ( - x27 + ( - x29 + ( - x31 + ( - x32))))))))))) * x4^2)))) * x1 + ((( (-2/1) * x6 + ( - x15 + ( - x16 + ( - x17 + ( - x18 + ( - x19 + ( - x24 + ( - x27 + ( - x29 + ( - x31 + ( - x32))))))))))) * x2 + ((( - x6 + ( - x7 + ( - x8 + ( - x25 + ( - x26 + ( - x27 + ( - x29 + ( - x31 + ( - x32))))))))) * x4 + ( (2/1) * x6 + ( (2/1) * x7 + (x16 + (x17 + (x18 + (x19 + (x20 + (x21 + (x22 + (x23 + ( (2/1) * x24 + ( (2/1) * x27 + ( (2/1) * x29 + ( (2/1) * x31 + ( (2/1) * x32)))))))))))))))) * x3 + ((x6 + (x8 + (x17 + (x18 + (x19 + (x24 + (x27 + (x29 + (x31 + (x32)))))))))) * x4))) * x2 + ((( (-2/1) * x7 + ( - x21 + ( - x22 + ( - x23 + ( - x24 + ( - x27 + ( - x29 + ( - x31 + ( - x32))))))))) * x3 + ((x7 + (x8 + (x22 + (x23 + (x24 + (x27 + (x29 + (x31 + (x32))))))))) * x4)) * x3 + (( - x8 + ( - x32)) * x4)));

interval=[[4 6.36];[4 6.36];[4 6.36];[4 6.36];repmat([-1 1],nparam,1)];
qsdp = box_norm(q,vars,interval);
%% Computation of the power matrix and the coefficient list
% qsdp = eval(q_norm);
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


% %% I (&J) complex sparcisity pattern
pattern = {[1,2,4];[1,3,4];[2,3,4]};
P = nvar:(nvar+nparam);

Crep = cell(size(pattern,1)*size(P,2),1);
for i=1:size(P,2)
    for j=1:size(pattern,1)
        Crep{(i-1)*size(pattern,1) +j} = [pattern{j} P(i)];
    end
end


%% I (complex sparsisity pattern)
dlmwrite([path '_ic.dat'],'');
for i=1:size(Crep)
    dlmwrite([path '_ic.dat']',Crep{i},'-append');
end
%% J = I

%% Classical pattern
dlmwrite([path '_i.dat'],I);
dlmwrite([path '_j.dat'],J);


