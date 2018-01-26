

%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 6;
nparam = 21;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);
q = (( (-2/1) * x7 + ( - x20 + ( - x21 + ( - x22 + ( - x23 + ( - x24 + ( - x25 + ( - x26 + ( - x27))))))))) * x1 + ((x7 + (x8 + (x21 + (x22 + (x23 + (x24 + (x25 + (x26 + (x27))))))))) * x2 + ((x7 + (x9 + (x22 + (x23 + (x24 + (x25 + (x26 + (x27)))))))) * x3 + (( - x7 + ( - x10 + ( - x23 + ( - x24 + ( - x25 + ( - x26 + ( - x27))))))) * x4 + ((x7 + (x11 + (x24 + (x25 + (x26 + (x27)))))) * x5 + ((x7 + (x12 + (x25 + (x26 + (x27))))) * x6)))))) * x1 + ((( - x8 + ( - x9 + ( - x16 + ( - x17 + ( - x19 + ( - x27)))))) * x3 + ((x8 + (x11 + (x13 + (x15 + (x17 + (x19 + (x27))))))) * x5)) * x2 + (((x9 + (x12 + (x14 + (x15 + (x17 + (x19 + (x27))))))) * x6) * x3 + ((( - x11 + ( - x12 + ( - x18 + ( - x19 + ( - x27))))) * x6) * x5)));
interval=[[4 6.36];[4 6.36];[4 6.36];[4 6.36];[4 6.36];[4 6.36];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1];[-1 1]];
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

mkdir('kepler0real');
dlmwrite('kepler0real/kepler0real_g.dat',G);
dlmwrite('kepler0real/kepler0real_s.dat',system_info);
dlmwrite('kepler0real/kepler0real_c.dat',c);
dlmwrite('kepler0real/kepler0real_p.dat',p);


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
dlmwrite('kepler0real/kepler0real_i.dat',I);
dlmwrite('kepler0real/kepler0real_j.dat',J);


