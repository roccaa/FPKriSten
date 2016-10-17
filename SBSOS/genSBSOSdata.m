    %%*********************************************************
    %% This code generates the Sparse Bounded Sums Of Squares 
    %% hierarchy for SDPT3
    %%*********************************************************
    %% F: objective polynomial
    %% G: a cell array containing the constraint polynomials
    %% I: a cell array explaining sparsity in variables
    %% J: a cell array explaining sparsity in constraints
    %% dd: the degree of Lagrangian dual relaxation
    %% kk: the degree of SOS relaxation
    %%*********************************************************


    function [bblk,AAt,CC,bb,constant,time,recy,degfg,oA,ob,oC,nlambda,oA_op,ob_op] = genSBSOSdata(F,G,I,J,dd,kk)
    
    time=cell(4,2);
    time{1,1} ='Init, analysis, writing';
    time{2,1} ='Compare coefficients (sparse representation)';
    time{3,1} ='generate Hab';
    time{4,1} ='Compare coefficients (certificates)';
    time{5,1} ='Remove dependent constraints and fixed variables';
    time{6,1} ='Total time';

    tic;

    %% OPTIONS
    
    RemoveConstant = true;
    RemoveDependent = false;
    RemoveFixvar = true;
    RecoverData = true;

    %% Analyse problem data
    
    fprintf('\n This is SBSOS \n')
    
    np1 = size(F,2);
    n = np1-1;                                % # variables total
    degf = max(sum(F(:,1:n),2));              % deg objective function
    m = length(G);                            % # constraints total
    if m
        gdeg = zeros(m,1);                    % # degree of each constraint
        for j=1:m
            gpower = G{j}(:,1:n);
            gdeg(j) = max(sum(gpower,2));
        end
    else
        gdeg = 0;
    end

    maxdeg = max([degf,dd*max(gdeg),2*kk]);  % maximal degree total
    degfg = max(degf,max(gdeg));
    %% ADDED BY ALEXANDRE ROCCA TO ENSURE AN EVEN DEGREE FOR HANDELMANN ?
%     if mod(degfg,2) == 1
%      degfg = degfg+1 % <======= LEADS TO SOME BUGS ....
%     end
    %%
    degFl = maxdeg;

    p = length(I);                           % number of sparse blocks
    
    fprintf('\n n = %1.0f, degf = %1.0f, m = %1.0f, d = %1.0f, k = %1.0f, maxdeg = %1.0f ',n,full(degf),m,dd,kk,full(maxdeg))

    %% Remove constant term before generating problem

    constant = 0;

    if (RemoveConstant)
        indConstants = find(ismember(F(:,1:n),zeros(1,n),'rows'));
        if sum(indConstants) > 0
            constant = sum(F(indConstants,np1));
            F(indConstants,:) = [];
        end
    end


    %% Preparing constraint generation
    
    nI = cellfun('length',I);          % # variables in each block
    mJ = cellfun('length',J);          % # constraints in each block
    fprintf('\n max_clique = %1.0f ',max(nI))
    
    sDim = zeros(p,1);        % size psd matrix for each block
    svDim = zeros(p,1);       % size of psd matrix represented in svec format
    HABcell = cell(p,1);       % all h(alpha,beta) for each block
    lDim = zeros(p,1);         % # (alpha,beta) in each block
    flMons = [];              % list of monomials for sparse representation
    fDim = zeros(p,1);        % number of monomials in sparse representation in each block
    tDim = 1;

    time{1,2} = toc;
    tic
    fprintf('\n Generating Hab polynomials for %d certificates ',p)
    for j=1:p
        sDim(j) = nchoosek(nI(j) + kk, nI(j)); 
        svDim(j) = sDim(j)*(sDim(j)+1)/2;
        if m
            GG = cell(mJ(j),1);                % constraints active in this block
            for i=1:mJ(j)
                GG{i} = G{J{j}(i)};
            end
        [HABcell{j},lDim(j)] = hpols(GG,dd);
        end    
        newmons = SparseMons(n,degFl,I{j});
        flMons = [flMons; newmons]; %#ok
        fDim(j) = size(newmons,1);
    end

    time{3,2} = toc;
    tic

    fdim = sum(fDim);      % total number of sparse representation monomials
    ldim = sum(lDim);       % number of lagrange multipliers
    %% POSSIBLE MODIFICATION ALEXANDRE ROCCA
    if(kk==0)
         sDim = zeros(p,1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    sdim = sum(sDim);         % size psd blk
    svdim = sum(svDim);
    %% MODIF ROCCA
    svdim = 0; 
    sdim = 0;
    %%%%%
    tdim = 1 ;                % number of optimizing variables

    nvar = svdim + ldim + tdim + fdim; % total number of variables

    %% Generate Variables for sBSOS and Contributution to the optimization criterion

    if m
        blk = cell(4,2);
        C = cell(4,1);
        blk{2,1} = 'l';          % lagrange multipliers
        blk{2,2} = ldim;
        C{2,1} = zeros(ldim,1);  % do not contribute to the optimization criterion
        ct = 2;
    else
        blk = cell(3,2);         % no constraints, no lagrange multipliers
        C = cell(3,1);
        ct = 1;
    end

    blk{1,1} = 's';              % psd variables
    blk{1,2} = sDim;
    C{1,1} = zeros(sdim);        % do not contribute to the optimization criterion

    ct = ct+1;
    blk{ct,1} = 'u';              % optimization variables
    blk{ct,2} = 1;
    C{ct,1} = -1;                 % opitmization criterion is max t

    ct = ct+1;
    blk{ct,1} = 'u';              % variables for sparse representation
    blk{ct,2} = fdim;
    C{ct,1} = zeros(fdim,1);      % do not contribute to the optimization criterion

    time{1,2} = time{1,2}+ toc;
    tic

    %% Generate Coefficients for Constraints

    %% Sparsity
    fprintf('\n Comparing coefficients of representing polynomials ...  ')
    Fexp = F(:,1:n);
    Fcof = F(:,np1);
    Flcof = ones(fdim,1);
    var = [spdiags([1;Flcof],0,fdim+1,fdim+2); sparse(size(F,1),fdim+1),Fcof];
    var_op = [spdiags([1;Flcof],0,fdim+1,fdim+2); sparse(size(F,1),fdim+1),-Fcof]; %% Added Alexandre Rocca
    
    
    tmp = compCoef([sparse(1,n);flMons;Fexp],var);
    sparsityb = tmp(:,fdim+2);
    sparsityAt = tmp(:,1:fdim+1); 
    
    %% ADDED ROCCA
    tmp_op  = compCoef([sparse(1,n);flMons;Fexp],var_op);
    sparsityb_op = tmp_op(:,fdim+2);
    sparsityAt_op = tmp_op(:,1:fdim+1); 
    %%
    
    time{2,2} = toc;
    tic

    %% Handelmann
    
    fprintf('\n Comparing coefficients for %d certificates ... ',p)
    sinit = 0;
    linit = svdim;
    finit = svdim + ldim + tdim;
    flinit = 0;
    certAt = [];
    certAt_op = [];
    %%%%%%%%%%%%%%%%%%%%%%
    for i= 1:p
        if m
            ind =cell2mat(cellfun(@size,HABcell{i},'uni',0));
            ind = ind(:,1);
            nhexp = sum(ind);
        else
            nhexp = 0;
        end

        var = sparse(svDim(i) + nhexp + tDim +fDim(i),nvar);
        %% ADDED ALEXANDRE ROCCA
        var_op =  sparse(svDim(i) + nhexp + tDim +fDim(i),nvar);
        
        %% ANCIENNE VERSION
%         [sexp, scof] = sospol(I{i},n,kk);
        %% AJOUT ALEXANDRE ROCCA
        % Eliminating the computation of the sos polynomial if k = 0
        if(kk~=0)
            [sexp, scof] = sospol(I{i},n,kk);
        else
            sexp = zeros(1,n);
            scof = 0;
        end
        %% ############
        for j = 1:svDim(i)
            var(j,sinit + j) = -scof(j);%#ok
            var_op(j,sinit + j) = -scof(j); %% ALEXANDRE ROCCA
        end
        mct = svDim(i);
        if m
            tmp = cell2mat(HABcell{i});
            habexp = tmp(:,1:n);
            habcof = tmp(:,np1);
            cct = 0;
            %% HERE eliminate loop ? (ROCCA)
    
            %% New vectorized version
            temp_C = cell2mat(arrayfun(@repmat, 1+linit:linit+lDim(i), ones(1,size(ind,1)), ind','UniformOutput',0));
            [rows,cols,val] = find(var);   
            var = sparse([rows mct+1:mct+sum(ind)],[cols temp_C],[val -habcof(cct+1:cct+sum(ind))],size(var,1),size(var,2));
            %% ADDED for min/max
            [rows,cols,val] = find(var_op);   
            var_op = sparse([rows mct+1:mct+sum(ind)],[cols temp_C],[val -habcof(cct+1:cct+sum(ind))],size(var_op,1),size(var_op,2));
            %%%%%
            mct = mct+sum(ind);
            cct = cct+sum(ind);
  
            %% OLD VERSION 
%             for j = 1:lDim(i)
%                 for k = 1:ind(j)
%                     mct = mct +1;
%                     cct = cct+1;
%                     var(mct,j+linit) = -habcof(cct);
%                 end
%             end
         
%            assert(prod(isequal(var_copy,var))>0)

        end
        
        fexp = flMons(flinit+1:flinit+fDim(i),:);
        flinit = flinit+fDim(i);
        
        for j =1:fDim(i)
            mct = mct +1;
            var(mct, finit + j) = 1;
            var_op(mct, finit + j) = 1;
        end
        if m
            exponents = [sexp;habexp;fexp];
        else
            exponents = [sexp;fexp];
        end
        
        
        at = compCoef(exponents,var);
        at_op = compCoef(exponents,var_op);  %% ROCCA
%         pause
%         size(at)
%         pause
        certAt = [certAt; at];
        certAt_op =  [certAt_op; at_op]; %% ROCCA
        
        sinit = sinit + svDim(i);
        linit = linit + lDim(i);
        finit = finit + fDim(i);
    end

    certb = zeros(size(certAt,1),1);
    time{4,2} = toc;
    tic

    %% Write constraints
    
    certAt = certAt';
    certAt_op = certAt_op'; %% ROCCA
    sparsityAt = sparsityAt';
    sparsityAt_op = sparsityAt_op'; %% ROCCA
    nspcons = size(sparsityAt,2);
    sparsityAt = [sparse(svdim,nspcons);sparse(ldim,nspcons);sparsityAt];
    sparsityAt_op = [sparse(svdim,nspcons);sparse(ldim,nspcons);sparsityAt_op]; %% ROCCA
    nspcons = nspcons - 1;
    
    allAt = [certAt,sparsityAt];  
    allAt_op = [certAt_op,sparsityAt_op];  %% ROCCA
    b = [certb;sparsityb];
    b_op = [certb;sparsityb_op]; %% ROCCA
%     pause
%     size(certb)
%     size(sparsityb)
%     pause

    if m
        At=cell(4,1);
        At_op = cell(4,1); %% ROCCA
    else
        At=cell(3,1);
    end

    At{1} = allAt(1:svdim,:);
    At_op{1} = allAt_op(1:svdim,:); %% ROCCA
    allAt(1:svdim,:) = [];
    allAt_op(1:svdim,:) = [];%% ROCCA
    if m 
    At{2} = allAt(1:ldim,:);
    allAt(1:ldim,:) = [];
    At_op{2} = allAt_op(1:ldim,:); %% ROCCA
    allAt_op(1:ldim,:) = []; %% ROCCA
    ct = 3;
    else
        ct =2;
    end
    
    At{ct} = allAt(1:tdim,:);
    allAt(1:tdim,:) = [];
    At{ct+1} = allAt;

    At_op{ct} = allAt_op(1:tdim,:); %% ROCCA
    allAt_op(1:tdim,:) = []; %% ROCCA
    At_op{ct+1} = allAt_op; %% ROCCA


    time{1,2} = time{1,2} + toc;
%     pause
%     disp('blk');
%     size(blk)
%     blk
%     disp('At');
%     size(At)
%     At
%     disp('C');
%     size(C)
%     C
%     disp('b');
%     size(b)
%     pause



  
    %% Recuparating info for LP only
    if(kk==0)
        oA = cell2mat(At(2:size(At,1)))';
        oA_op = cell2mat(At_op(2:size(At_op,1)))';
        nlambda = size(At{2},1);
        oC = cell2mat(C(2:size(At,1)))';
        ob = b;
        ob_op = b_op;
%         time{5,2} = toc;
%         %% 
%         time{6,2} = time{1,2} +time{2,2} + time{3,2} + time{4,2} + time{5,2};
    end


    tic

    %% optimize SDP fo sdpt3
    fprintf('\n Optimize SDP for solver ... ')
    recy.I = I;
    recy.At = At;
    recy.degFl = degFl;
    recy.nspcons = nspcons;
    recy.RemoveFixvar = RemoveFixvar;
    recy.RemoveDependent = RemoveDependent;
    
    % remove fixed variables

    if RemoveFixvar
        [blk1,At1,C1,b1,R,scols,m] = removefixvar(blk,At,C,b,nspcons,RecoverData);
        recy.R = R; recy.scols = scols; recy.m = m; 
    else
        blk1 = blk;
        At1 = At;
        C1 = C;
        b1 = b;
    end
       
    % remove dependent constraints    
    if RemoveDependent
        y0 = zeros(length(b1),1);
        [AAt,bb,y0,idxB,neardepconstr,feasible] =  checkdepconstr(blk1,At1,b1,y0,1);
        fprintf('\n Removed %2.0f nearly dependent constraints',neardepconstr)
        CC = C1;
        bblk = blk1;
        recy.y0 = y0; recy.idxB = idxB;
    else
        bblk = blk1;
        AAt = At1;
        CC = C1;
        bb = b1;
    end
    
   

    time{5,2} = toc;
    
    %% 
    time{6,2} = time{1,2} +time{2,2} + time{3,2} + time{4,2} + time{5,2};
%     time
    %*************************************************************************%
    end