
%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
complex_sparse = 0;
nvar = 1;
name = 'simple_test';
nparam = 3;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);

% q = (x1+1)*x1*(x2 +x3);
q = (x1*x1-x1)*(x4)+x1*x1*x3+(2*x1*x1-x1)*x2;

interval=[[0 1];repmat([-1 1],nparam,1)];
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

path = [name '/' name];
mkdir(name);
dlmwrite([path '_g.dat'],G);
dlmwrite([path '_s.dat'],system_info);
dlmwrite([path '_c.dat'],c);
dlmwrite([path '_p.dat'],p);


% %% I (&J) complex sparcisity pattern
% pattern = {[1,2,4];[1,3,4];[2,3,4]};
% P = nvar:(nvar+nparam);
% 
% Crep = cell(size(pattern,1)*size(P,2),1);
% for i=1:size(P,2)
%     for j=1:size(pattern,1)
%         Crep{(i-1)*size(pattern,1) +j} = [pattern{j} P(i)];
%     end
% end
% 
% 
% %% I (complex sparsisity pattern)
% dlmwrite([path '_ic.dat'],'');
% for i=1:size(Crep)
%     dlmwrite([path '_ic.dat']',Crep{i},'-append');
% end
%% J = I

%% Classical pattern
dlmwrite([path '_i.dat'],I);
dlmwrite([path '_j.dat'],J);


