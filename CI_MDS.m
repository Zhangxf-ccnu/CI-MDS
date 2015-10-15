function [X,gamma_G] = CI_MDS(PPI_file, ell, output_file, method)

% CI_MDS is the corresponding matlab function of the proposed model. It reads data (PPI network) from text
% file (PPI_file) using function: PPI_read, and determines driver proteins using the model
% described in the paper. It also writes the identified driver proteins
% into file 'output_file' where each line is a identified driver protein.


% Inputs:
%   PPI_file: the input file name of the PPI network, where each line
%   contains two proteins defining the interaction between them.
%   Example: ABI1	ABL1

%   ell£ºthe distant parameter used for calculating collect influence of proteins. The default value is 1.

%   output_file: the file into which the identified driver proteins will be written. Each line is a
%   predicted driver protein. The default value is 'output_file.txt'.
%   Example: ABL1

%   method: the optimization methods ('lp_solve' or 'intlinprog') that will be used to solve the binary
%   integer-programming problem. The default value is 'intlinprog'.

% Outputs:
%   X: the binary vector that indicates whether a protein is a member of
%   detetermined MDS.

%   gamma_G: the domination number

if nargin < 4
    method = 'intlinprog';
end

if nargin < 3
    output_file = 'output_file.txt';
end

if nargin < 2
    ell = 1;
end


%read PPI data 
Data_set = PPI_read(PPI_file);
A = Data_set.PPI;

fprintf('Computing collective influence ...')
fprintf('\n')
%compute the number of proteins
n = size(A,1);
%compute collective influence of proteins
A = A - diag(diag(A));
Degree = sum(A,2);
shortest = graphallshortestpaths(A);
Ball = shortest == ell;
% EQ. (3)
CI  = (Degree - 1).*(Ball*(Degree -1 ));

% add self-interactions
A = A-diag(diag(A)) +eye(size(A));

% solve the standard MDS model EQ.(1)
fprintf('Solve standard MDS model ...')
fprintf('\n')

switch method
    % choose optimization method 'lp_solve'
    case 'lp_solve'
        [~,X,~,~] = lp_solve(-ones(n,1) ,A, ones(n,1),ones(n,1),zeros(n,1),ones(n,1),1:n);
    case 'intlinprog'
        X = intlinprog(ones(n,1), 1:n, -A, -ones(n,1), [], [], zeros(n,1), ones(n,1));
end


% computing the domination number accoring to EQ. (2)
fprintf('Computing domination number ...')
fprintf('\n')
gamma_G = sum(X);

% solve the  CI-MDS model EQ.(4)
fprintf('Slove the CI-MDS model ...')
fprintf('\n')

switch method
    % choose optimization method 'lp_solve'
    case 'lp_solve'
        [~,X,~,~] = lp_solve( CI, [A; ones(1,n)], [ones(n,1); gamma_G],[ones(n,1);0],zeros(n,1),ones(n,1),1:n);
    case 'intlinprog'
        X = intlinprog( -CI, 1:n, -A, -ones(n,1), ones(1,n), gamma_G, zeros(n,1), ones(n,1));
end

% Write the result of predicted driver proteins into file 'output_file'.
fprintf('Print the results ...')
fprintf('\n')
result_print(output_file, Data_set.Protein,  X)



function Data_set = PPI_read(PPI_profie) %#ok<*DEFNU>
fprintf('Reading data...')
fprintf('\n')
fid_ppi=fopen(PPI_profie);
temp_PPI=textscan(fid_ppi,'%s%s%*[^\n]','delimiter','\t');
fclose(fid_ppi);

Data_set.Protein = union(temp_PPI{1},temp_PPI{2});
[~,Locb_1] = ismember(temp_PPI{1}, Data_set.Protein);
[~,Locb_2] = ismember(temp_PPI{2}, Data_set.Protein);
Data_set.PPI = sparse(Locb_1,Locb_2,1,length(Data_set.Protein),length(Data_set.Protein));
Data_set.PPI = Data_set.PPI + Data_set.PPI';
Data_set.PPI = Data_set.PPI- diag(diag(Data_set.PPI));


function result_print(result_file, Proteins,  X)
Driver_proteins = Proteins(logical(X));
fid = fopen(result_file,'w');
for i = 1:length(Driver_proteins)
    fprintf(fid, '%s',Driver_proteins{i});
    fprintf(fid, '\n');
end
fclose(fid);
fprintf(['The predicted driver proteins have been written into file ', result_file])
fprintf('\n')
