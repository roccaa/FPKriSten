%% Solve a given example
%%
%% Output ==> ?
function [ bound,build_time,solving_time ] = solve_examples( F,G,I,J,d,k,tag )

[blk,At,C,b,constant,time,recy,degfg,oA,ob,oC,nlambda,oA_op,ob_op] = genSBSOSdata(F,G,I,J,d,k);

% blk
% size(blk)

% [obj,X,y,Z,info] = sdpt3(blk,At,C,b);
% if(info.termcode~=0)
%     [obj,X,y,Z,info] = sqlp(blk,At,C,b);
% end
% 
% if(strcmp(tag,'MAX'))
%     bound = -(-obj(1)+constant) 
% else
%     bound = (-obj(1)+constant) 
% end
% 

% 
% 
% disp('VERIFICATION OF THE RANK');
% YY = recovery(recy,y);
% [RANK,mineig] = rankSBSOS(YY,I,degfg)
% rnk = sum(RANK)/length(RANK)

build_time = time{6,2};
time

% disp('Linprog solving');

% tic
% [x,fval] = linprog(oC,[],[],oA,ob,sparse(nlambda,1));
% toc
% fval
% if(strcmp(tag,'MAX'))
%     bound = -(-fval+constant) 
% else
%     bound = (-fval+constant) 
% end
% 
% disp('GLPKMEX solving');
% 
% ctype = repmat('S',1,size(oA,1)); % All the constraint are equality constraints
% vartype = repmat('C',1,size(oA,2)); % All the varaible are 'continuous' variables
% 
% s = 1; % minimization pb
% 
% param.msglev = 1; % All messages
%  param.itlim = 100; % max 100 iterations
% % param.lpsolver = 1; % (1) for simplex, (2) for interior point
% 
% FLOAT_MIN  = -realmax('single'); % Borne min = -MAX_FLOAT
% FLOAT_MIN = double(FLOAT_MIN);
%  lb = [sparse(nlambda,1);repmat(FLOAT_MIN,size(oA,2)-nlambda,1)];
%  % non-negative lambda
% 
% tic
%   [xopt, fmin, status, extra] = glpk(full(oC), oA, full(ob), full(lb), [], ctype, vartype,s,param);
% toc
% fmin
% status
% extra
% if(strcmp(tag,'MAX'))
%     bound = -(-fmin+constant) 
% else
%     bound = (-fmin+constant) 
% endFLOAT_MIN  = -realmax

oC = sparse(oC);
oA = sparse(oA);
% pause
% At
% size(oA)
% blk
% size(b)
% pause
% ob = sparse(ob);
ob_op = sparse(ob_op);
oA_op = sparse(oA_op);
FLOAT_MIN  = -realmax
lb = [sparse(nlambda,1);repmat(FLOAT_MIN,size(oA,2)-nlambda,1)];

options = cplexoptimset('Algorithm','interior-point','Display','off');

disp('CPLEX solving (MAX)');
tic 
[x fval exitflag] = cplexlp(oC,[],[],oA,ob,lb,[],[],options);
ts1 = toc
exitflag
fval;
bound1 = -(-fval+constant)

disp('CPLEX solving (MIN)');
tic 
[x fval exitflag] = cplexlp(oC,[],[],oA_op,ob_op,lb,[],[],options);
ts2 = toc
exitflag
total_time  = build_time + ts1 + ts2
solving_time = ts1 +ts2;
fval;
bound2 = (-fval-constant)
bound = max(abs(bound1),abs(bound2))


% 
% pause
% oC
% pause
% oA
% pause
% ob
% pause
% lb
% pause
% ctype
% pause
% vartype


end

