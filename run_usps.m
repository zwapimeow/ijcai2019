%% Settings
% setup
nCluster = 5; % number of subjects.

nSample = 250; % number of images for each digit

% normalization
normalizeColumn = @(data) cnormalize_inplace(data);
% representation     
% buildRepresentation = @(data) ssc_relaxed(data, 0.1);
% buildRepresentation = @(data) SSCOMP(data, 8, 1e-6); % second parameter is sparsity
buildRepresentation = @(data) oursSSCOMP(data, 8, 1e-6); % second parameter is sparsity
% buildRepresentation = @(data) SSCROMP(data, 8, 1e-6); % second parameter is sparsity
% buildRepresentation = @(data) oursSSCROMP(data, 8, 1e-6); % second parameter is sparsity
% spectral clustering       
genLabel = @(affinity, nCluster) SpectralClustering(affinity, nCluster, 'Eig_Solver', 'eigs');
% genLabel = @(affinity, nCluster) ncutW(affinity, nCluster);

%% Load data
% Data is preprocessed and saved in the .mat file. 
% usps_data.mat : size = [(16*16),(10*1100)] 
% use 'im2double' make data from uint8 to double
% imagesize 16*16 , 10 numbers(1-9;0) * 1100 images/number = 1100 
% usps_lable.mat : size = [(10*1100),1] , (1-9;0)*1100
load usps_data.mat 
load usps_label.mat
N_subject = 10;

%% Clustering
nExperiment = 20;
results = zeros(nExperiment, 7); %results

nnzR = 0;
nnzA = 0;
meanK = 0;

for iExperiment = 1:nExperiment
%     digit_set = 0:9; % set of digits to test on, e.g. [2, 0]. Pick randomly if empty.   
    digit_set = [];
    % prepare data
    if isempty(digit_set)
        rng(iExperiment); Digits = randperm(10, nCluster) - 1;
    else
        Digits = digit_set;
    end
    if length(nSample) == 1
        nSample = ones(1, nCluster) * nSample;
    end
    mask = zeros(1, sum(nSample));
    s = zeros(1, sum(nSample));
    nSample_cum = [0, cumsum(nSample)];
    for iK = 1:nCluster % randomly take data for each digit.
        allpos = find( usps_label == Digits(iK) );
        rng( (iExperiment-1) * 10 + iK );
        selpos = allpos( randperm(length(allpos), nSample(iK)) );

        mask( nSample_cum(iK) + 1 : nSample_cum(iK+1) ) = selpos;
        s( nSample_cum(iK) + 1 : nSample_cum(iK+1) ) = iK * ones(1, nSample(iK));
    end
    X = usps_data(:, mask);
    N = length(s);
    
%     X = imnoise(X,'salt',0.4);
        
    % Clustering
    tic;
    
    % normalization
%     fprintf('Normalization...\n')
    X = normalizeColumn(X);
    % generate representation
%     fprintf('Representation...\n')

%     X = awgn(X, 10*log10(128) );

%     R = buildRepresentation(X);
%     K=[];
    [R,K] = buildRepresentation(X);
    time1 = toc;
    % generate affinity
%     fprintf('Affinity...\n')
    R(1:N+1:end) = 0;
    R = cnormalize(R, Inf);
    A = abs(R) + abs(R)';
    % generate label
%     fprintf('Generate label...\n')
    groups = genLabel(A, nCluster);   
%     [~,groups]=max(groups,[],2);% if genLable~=ncutW then % this sentence

    time2 = toc;
    
    % Evaluation
    perc = evalSSR_perc( R, s );
    ssr = evalSSR_error( R, s );   
    conn = evalConn( A, s);
    accr  = evalAccuracy(s, groups);
    nmi_per = nmi(s,groups);

    % output
    dataformat = '%d-th experiment: perc = %f, ssr = %f, conn = %f, accr = %f, time1 = %f, time2 = %f, nmi = %f\n';
    dataValue = [iExperiment, perc, ssr, conn, accr, time1, time2-time1];
    fprintf(dataformat, dataValue, nmi_per);
    % record
    results(iExperiment, :) = dataValue;
    nnzR = nnzR + nnz(R);
    nnzA = nnzA + nnz(A);
    meanK = meanK + mean(K);
end
% output
dataValue = mean(results, 1);
fprintf('\nAverage: perc = %f, ssr = %f, conn = %f, accr = %f, time1 = %f, time2 = %f, nmi = %f\n', dataValue(2:end), nmi_per);
fprintf("%f , %f , %f\n", nnzR/nExperiment, nnzA/nExperiment, meanK/nExperiment);
