
%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
use_semialgebraic = 1;
complex_sparse = 0;
nvar = 6;
name ='floudas3_3';
nparam = 25;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);

% Model
q = (( -50 * x7 + ( -25 * x13 + ( -50 * x14 + ( -25 * x15 + ( -25 * x16 + ( -25 * x19 + ( -25 * x22 + ( -25 * x25 + ( -25 * x28 + ( -25 * x31)))))))))) * x1 + ( 100 * x7 + ( 100 * x13 + ( 200 * x14 + ( 100 * x15 + ( 100 * x16 + ( 100 * x19 + ( 100 * x22 + ( 100 * x25 + ( 100 * x28 + ( 100 * x31))))))))))) * x1 + ((( -2 * x8 + ( -2 * x17 + ( - x18 + ( - x19 + ( - x22 + ( - x25 + ( - x28 + ( - x31)))))))) * x2 + ( 4 * x8 + ( 8 * x17 + ( 4 * x18 + ( 4 * x19 + ( 4 * x22 + ( 4 * x25 + ( 4 * x28 + ( 4 * x31))))))))) * x2 + ((( -2 * x9 + ( -2 * x20 + ( - x21 + ( - x22 + ( - x25 + ( - x28 + ( - x31))))))) * x3 + ( 2 * x9 + ( 4 * x20 + ( 2 * x21 + ( 2 * x22 + ( 2 * x25 + ( 2 * x28 + ( 2 * x31)))))))) * x3 + ((( -2 * x10 + ( -2 * x23 + ( - x24 + ( - x25 + ( - x28 + ( - x31)))))) * x4 + ( 8 * x10 + ( 16 * x23 + ( 8 * x24 + ( 8 * x25 + ( 8 * x28 + ( 8 * x31))))))) * x4 + ((( -2 * x11 + ( -2 * x26 + ( - x27 + ( - x28 + ( - x31))))) * x5 + ( 2 * x11 + ( 4 * x26 + ( 2 * x27 + ( 2 * x28 + ( 2 * x31)))))) * x5 + ((( -2 * x12 + ( -2 * x29 + ( - x30 + ( - x31)))) * x6 + ( 8 * x12 + ( 16 * x29 + ( 8 * x30 + ( 8 * x31))))) * x6 + ( -100 * x13 + ( -200 * x14 + ( -100 * x15 + ( -100 * x16 + ( -8 * x17 + ( -4 * x18 + ( -104 * x19 + ( -2 * x20 + ( - x21 + ( -105 * x22 + ( -32 * x23 + ( -16 * x24 + ( -121 * x25 + ( -2 * x26 + ( - x27 + ( -122 * x28 + ( -32 * x29 + ( -16 * x30 + ( -138 * x31))))))))))))))))))))))));
cstr = [((x3-3)^2+x4-4)/6;
        ((x5-3)^2+x6-4)/10;
        (2-x1+3*x2)/20;
        (2+x1-x2)/8;
        (6-x1-x2)/6;
        (x1+x2-2)/10];

n_semialg =  size(cstr,1);
    
interval=[[0 6];[0 6];[1 5];[0 6];[1 5];[0 10];repmat([-1 1],nparam,1)];
qsdp = box_norm(q,vars,interval);

% normalizing semi-algebraic equation
for i=1:n_semialg
    cstr(i) = box_norm(cstr(i),vars,interval);
end
%% Computation of the power matrix and the coefficient list
% qsdp = eval(q_norm);
[powers,coefficients] = getexponentbase(qsdp,vars);
p = str2double(sdisplay(powers));
c = str2double(sdisplay(coefficients));



n = nvar+nparam;
[I,J] = build_box_sparcity(nvar,nparam);
% Completing J with semi-algebraic constraints
ci = nvar+nparam+1:nvar+nparam+n_semialg;

C = repmat(ci,nparam,1);
J = [J C];

system_info = [n complex_sparse n_semialg];
G = create_unitBox(n);

%% Similar to box only
path = [name '/' name];
mkdir(name);
dlmwrite([path '_s.dat'],system_info);
dlmwrite([path '_c.dat'],c);
dlmwrite([path '_p.dat'],p);
%% Classical pattern
dlmwrite([path '_i.dat'],I);
dlmwrite([path '_j.dat'],J);


%% getting info on Semialgebraic
for i=1:n_semialg
    [powersloc,coefficientsloc] = getexponentbase(cstr(i),vars);
     ploc = str2double(sdisplay(powersloc));
     cloc = str2double(sdisplay(coefficientsloc));
     dlmwrite([path '_p_Gsa_' num2str(i) '.dat'],ploc);
     dlmwrite([path '_c_Gsa_' num2str(i) '.dat'],cloc);
end

dlmwrite([path '_g.dat'],G);


