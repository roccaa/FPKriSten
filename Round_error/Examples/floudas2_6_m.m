
%% Buildind the model files for kepler0 model
%% Normalization (from box representation) in [0:1]
%% Use the symbolic toolbox
use_semialgebraic = 1;
complex_sparse = 0;
nvar = 10;
name ='floudas2_6';
nparam = 50;
[str,vars]  = build_sdpvar(nvar,nparam);
eval(str);
vars = eval(vars);

% Model
q = (( -100 * x11 + ( -50 * x40 + ( -50 * x42 + ( -50 * x44 + ( -50 * x46 + ( -50 * x48 + ( -50 * x50 + ( -50 * x52 + ( -50 * x54 + ( -50 * x56 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60))))))))))))) * x1 + ( 48 * x11 + ( 48 * x21 + ( 48 * x23 + ( 48 * x25 + ( 48 * x27 + ( 48 * x29 + ( 48 * x31 + ( 48 * x33 + ( 48 * x35 + ( 48 * x37 + ( 48 * x39 + ( 48 * x60))))))))))))) * x1 + ((( -100 * x12 + ( -50 * x41 + ( -50 * x42 + ( -50 * x44 + ( -50 * x46 + ( -50 * x48 + ( -50 * x50 + ( -50 * x52 + ( -50 * x54 + ( -50 * x56 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60))))))))))))) * x2 + ( 42 * x12 + ( 42 * x22 + ( 42 * x23 + ( 42 * x25 + ( 42 * x27 + ( 42 * x29 + ( 42 * x31 + ( 42 * x33 + ( 42 * x35 + ( 42 * x37 + ( 42 * x39 + ( 42 * x60))))))))))))) * x2 + ((( -100 * x13 + ( -50 * x43 + ( -50 * x44 + ( -50 * x46 + ( -50 * x48 + ( -50 * x50 + ( -50 * x52 + ( -50 * x54 + ( -50 * x56 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60)))))))))))) * x3 + ( 48 * x13 + ( 48 * x24 + ( 48 * x25 + ( 48 * x27 + ( 48 * x29 + ( 48 * x31 + ( 48 * x33 + ( 48 * x35 + ( 48 * x37 + ( 48 * x39 + ( 48 * x60)))))))))))) * x3 + ((( -100 * x14 + ( -50 * x45 + ( -50 * x46 + ( -50 * x48 + ( -50 * x50 + ( -50 * x52 + ( -50 * x54 + ( -50 * x56 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60))))))))))) * x4 + ( 45 * x14 + ( 45 * x26 + ( 45 * x27 + ( 45 * x29 + ( 45 * x31 + ( 45 * x33 + ( 45 * x35 + ( 45 * x37 + ( 45 * x39 + ( 45 * x60))))))))))) * x4 + ((( -100 * x15 + ( -50 * x47 + ( -50 * x48 + ( -50 * x50 + ( -50 * x52 + ( -50 * x54 + ( -50 * x56 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60)))))))))) * x5 + ( 44 * x15 + ( 44 * x28 + ( 44 * x29 + ( 44 * x31 + ( 44 * x33 + ( 44 * x35 + ( 44 * x37 + ( 44 * x39 + ( 44 * x60)))))))))) * x5 + ((( -100 * x16 + ( -50 * x49 + ( -50 * x50 + ( -50 * x52 + ( -50 * x54 + ( -50 * x56 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60))))))))) * x6 + ( 41 * x16 + ( 41 * x30 + ( 41 * x31 + ( 41 * x33 + ( 41 * x35 + ( 41 * x37 + ( 41 * x39 + ( 41 * x60))))))))) * x6 + ((( -100 * x17 + ( -50 * x51 + ( -50 * x52 + ( -50 * x54 + ( -50 * x56 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60)))))))) * x7 + ( 47 * x17 + ( 47 * x32 + ( 47 * x33 + ( 47 * x35 + ( 47 * x37 + ( 47 * x39 + ( 47 * x60)))))))) * x7 + ((( -100 * x18 + ( -50 * x53 + ( -50 * x54 + ( -50 * x56 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60))))))) * x8 + ( 42 * x18 + ( 42 * x34 + ( 42 * x35 + ( 42 * x37 + ( 42 * x39 + ( 42 * x60))))))) * x8 + ((( -100 * x19 + ( -50 * x55 + ( -50 * x56 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60)))))) * x9 + ( 45 * x19 + ( 45 * x36 + ( 45 * x37 + ( 45 * x39 + ( 45 * x60)))))) * x9 + ((( -100 * x20 + ( -50 * x57 + ( -50 * x58 + ( -50 * x59 + ( -50 * x60))))) * x10 + ( 46 * x20 + ( 46 * x38 + ( 46 * x39 + ( 46 * x60))))) * x10)))))))));
cstr = [(-4+2*x1+6*x2+1*x3+0*x4+3*x5+3*x6+2*x7+6*x8+2*x9+2*x10)/23;
        (22-(6*x1-5*x2+8*x3-3*x4+0*x5+1*x6+3*x7+8*x8+9*x9-3*x10))/33;
        (-6-(5*x1+6*x2+5*x3+3*x4+8*x5-8*x6+9*x7+2*x8+0*x9-9*x10))/11;
        (-23-(9*x1+5*x2+0*x3-9*x4+1*x5-8*x6+3*x7-9*x8-9*x9-3*x10))/23;
        (-12-(-8*x1+7*x2-4*x3-5*x4-9*x5+1*x6-7*x7-1*x8+3*x9-2*x10))/24];

n_semialg =  size(cstr,1);
    
interval=[repmat([0 1],nvar,1);repmat([-1 1],nparam,1)];
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


