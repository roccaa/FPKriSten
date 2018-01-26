
%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 6;
nparam = 15;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q = (( - x14 + ( - x15 + ( - x16 + ( - x17 + ( - x18 + ( - x19 + ( - x20 + ( - x21)))))))) * x1 + ((x15 + (x16 + (x17 + (x18 + (x19 + (x20 + (x21))))))) * x2 + ((x16 + (x17 + (x18 + (x19 + (x20 + (x21)))))) * x3 + (( - x17 + ( - x18 + ( - x19 + ( - x20 + ( - x21))))) * x4 + ((x18 + (x19 + (x20 + (x21)))) * x5 + ((x19 + (x20 + (x21))) * x6)))))) * x1 + ((( - x10 + ( - x11 + ( - x13 + ( - x21)))) * x3 + ((x7 + (x9 + (x11 + (x13 + (x21))))) * x5)) * x2 + (((x8 + (x9 + (x11 + (x13 + (x21))))) * x6) * x3 + ((( - x12 + ( - x13 + ( - x21))) * x6) * x5)));
interval=[[4 6.36];[4 6.36];[4 6.36];[4 6.36];[4 6.36];[4 6.36];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1]];
qsdp = box_norm(q,vars,interval);
%% Computation of the power matrix and the coefficient list
% qsdp = eval(q_norm);
[powers,coefficients] = getexponentbase(qsdp,vars);
p = str2double(sdisplay(powers));
c = str2double(sdisplay(coefficients));

n = nvar+nparam;
[I,J] = build_box_sparcity(nvar,nparam);
cstr = []
n_semialg =  size(cstr,1);
system_info = [n complex_sparse n_semialg];
G = create_unitBox(n);

mkdir('kepler0');
dlmwrite('kepler0/kepler0_g.dat',G);
dlmwrite('kepler0/kepler0_s.dat',system_info);
dlmwrite('kepler0/kepler0_c.dat',c);
dlmwrite('kepler0/kepler0_p.dat',p);


% %% I (&J) complex sparcisity pattern
% C = {[1,4];[1,2,3];[1,2,5];[1,5,6];[1,3,6]};
% P = [7,8,9,10,11,12,13,14,15,16,17,18,19,20,21];
% 
% Crep = cell(size(C,1)*size(P,2),1);
% for i=1:size(P,2)
%     for j=1:size(C,1)
%         Crep{(i-1)*size(C,1) +j} = [C{j} P(i)];
%     end
% end




%% I (complex sparsisity pattern)
% dlmwrite('kepler0/kepler0_i.dat','');
% for i=1:size(Crep)
%     dlmwrite('kepler0/kepler0_i.dat',Crep{i},'-append');
% end
% %% J
% dlmwrite('kepler0/kepler0_j.dat','');
% for i=1:size(Crep)
%     dlmwrite('kepler0/kepler0_j.dat',Crep{i},'-append');
% end
dlmwrite('kepler0/kepler0_i.dat',I);
dlmwrite('kepler0/kepler0_j.dat',J);


