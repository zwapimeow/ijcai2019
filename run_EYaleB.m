% This code tests the performance of SSC-OMP on the EYaleB face database.
% The code generates results in Table 2 of the paper
% C. You, D. Robinson, R. Vidal, Scalable Sparse Subspace Clustering by 
% Orthogonal Matching Pursuit, CVPR 2016.

% In this code, we apply SSC-OMP to the face images of a randomly picked 
% n \in \{2, 10, 20, 30, 38} subjects in the Extended Yale B database. Each
% subject has 64 images under different illumination conditions. 

% Instructions for running the code:
% - Download code for computing clustering accuracy. Go to
% http://www.cad.zju.edu.cn/home/dengcai/Data/Clustering.html and download
% bestMap.m and Hungarian.m. Alternatively, you can use your own function
% by redefining the function evalAccuracy.m
% - Run. You can modify the parameter "nCluster" below to run for different
% number of subjects. 

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

%% Settings
% setup
nCluster = 5; % number of subjects. <= 38

% dimension reduction
% reduceDimension = @(data) dimReduction_PCA(data, 500);

% normalization
normalizeColumn = @(data) cnormalize_inplace(data);

% representation
% buildRepresentation = @(data) ssc_relaxed(data, 0.1);
% buildRepresentation = @(data) SSCOMP(data, 8, 1e-6); % second parameter is sparsity
% buildRepresentation = @(data) oursSSCOMP(data, 8, 1e-6); % second parameter is sparsity
% buildRepresentation = @(data) SSCROMP(data, 8, 1e-6); % second parameter is sparsity
% buildRepresentation = @(data) oursSSCROMP(data, 8, 1e-6); % second parameter is sparsity

% spectral clustering   
genLabel = @(affinity, nCluster) SpectralClustering(affinity, nCluster, 'Eig_Solver', 'eigs');

%% Load data
% Data is preprocessed and saved in the .mat file. 
% EYALEB_DATA is a D by N matrix. Each column is a face image and N = 38 subjects * 64 images/subject = 2414.
% Each image is downsampled from 192*168 to D = 48*42 = 2016.
% EYALEB_LABEL is a 1 by N vector. 
% Each entry is the label for the corresponding column in EYALEB_DATA.

load ExtendedYaleB.mat EYALEB_DATA EYALEB_LABEL
N_subject = length(unique(EYALEB_LABEL));
%% Clustering
nExperiment = 20;
results = zeros(nExperiment, 6); %results
% for r = 1.1:0.1:3
%     fprintf('\n r:%f\n',r);

nnzR = 0;
nnzA = 0;
meanK = 0;

for iExperiment = 1:nExperiment
    % prepare data
    rng(iExperiment * 38 + nCluster);
    subjectIdx = randperm(N_subject, nCluster); % select #nCluster subjects
    datapointIdx = find(ismember(EYALEB_LABEL, subjectIdx));
    X = double(EYALEB_DATA(:, datapointIdx)); % data
    s = EYALEB_LABEL(datapointIdx); % label
    N = length(s);
    % clustering
    
%     X = X/256;
%     X = imnoise(X,'gaussian',0.1);
    
    tic;
%     fprintf('Dimension reduction...\n')
%    [~,X] = reduceDimension(X);
    % normalization
%     fprintf('Normalization...\n')
    X = normalizeColumn(X);
    % generate representation
%     fprintf('Representation...\n')
%     R = buildRepresentation(X);
%     K=[];
    [R,K] = buildRepresentation(X);
    % generate affinity
%     fprintf('Affinity...\n')
    R(1:N+1:end) = 0;
%     % R = cnormalize(R, Inf);
%     R1 = ROMP_mat_thr_col_sec(R,1e-6);
%     R1(1:N+1:end) = 0;
%     A_0 = R*R1;
%     A_1 = abs(A_0) + abs(A_0)';
% %     A = abs(R) + abs(R)';
%     R1 = ROMP_mat_thr_col_sec(A,1e-6);
%     R1(1:N+1:end) = 0;
%     A_0 = abs(R1) + abs(R1)';
%     A_1 = A*A_0;
    A = abs(R) + abs(R)';%0
%     A = abs(corr(R));A = A+A';

    % generate label
%     fprintf('Generate label...\n')
    groups = genLabel(A, nCluster); %0 
%     [~,groups]=max(groups,[],2);% if genLable~=ncutW then % this sentence
%     groups = genLabel(A_1, nCluster);
    time = toc;
%     
    % Evaluation
    perc = evalSSR_perc( R, s );%0
    ssr = evalSSR_error( R, s );%0   
%    conn = evalConn_1( A, s);%0
    conn = evalConn( A, s);
    accr  = evalAccuracy(s, groups);%0
    nmi_per = nmi(s,groups);
% 
%     perc = evalSSR_perc( A_0, s );
%     ssr = evalSSR_error( A_0, s );   
%     conn = evalConn( A_1, s);
%     accr  = evalAccuracy(s, groups);

    % output
    dataformat = '%d-th experiment: perc = %f, ssr = %f, conn = %f, accr = %f, time = %f, nmi = %f\n';
    dataValue = [iExperiment, perc, ssr, conn, accr, time];
    fprintf(dataformat, dataValue, nmi_per);
%    nmi_per = nmi(s,groups);
%    fprintf('nmi: %f\n', nmi_per);
    % record
    results(iExperiment, :) = dataValue;
    nnzR = nnzR + nnz(R);
    nnzA = nnzA + nnz(A);
    meanK = meanK + mean(K);
end
% output
dataValue = mean(results, 1);
fprintf('\nAverage: perc = %f, ssr = %f, conn = %f, accr = %f, time = %f, nmi = %f\n', dataValue(2:end), nmi_per);
fprintf("%f , %f , %f\n", nnzR/nExperiment, nnzA/nExperiment, meanK/nExperiment);
fprintf("%f , %f , %f , %f\n", max(results(:,5)), mean(results(:,5)), min(results(:,5)), var(results(:,5)));
% end
