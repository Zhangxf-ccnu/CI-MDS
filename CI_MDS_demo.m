%%
clear 
% set ell = 1
ell = 1;


% CC-MDS model: optimization method "lp_solve"; ell =1; tissue = 'Adipose'
[X, domination_number] = CI_MDS('.\TissueNet\Adipose.txt', ell,'Adipose_MDS_lp_solve.txt', 'lp_solve');

% CC-MDS model: optimization method "intlinprog"; ell =1; tissue = 'Adipose'
[X, domination_number] = CI_MDS('.\TissueNet\Adipose.txt', ell,'Adipose_MDS_intlinprog.txt', 'intlinprog');

