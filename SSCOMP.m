function [C,Karray] = SSCOMP(X, K, thr)
    MEMORY_TOTAL = 0.1 * 10^9;
    [~, N] = size(X); 
    Xn = X;
    Xn = cnormalize(X);
    S = zeros(N, N);
    Val = zeros(N, K);
    res = Xn; 
    nonzero = zeros(N,1);

%     euDist = X'*X;                      %,,,,,,,
%     euDist = sort(euDist,2,'descend');  %,,,,,,,
%     euDist = euDist(:,2:K);             %,,,,,,,
%     euDist = mean(euDist,2);            %,,,,,,,
%     maxDist = max(euDist);              %,,,,,,,
%     minDist = min(euDist(euDist~=0));   %,,,,,,,
%     dist = maxDist - minDist;           %,,,,,,,
%     euDist = K*(euDist - minDist)/dist; %,,,,,,,
%     move = K - round(mean(euDist));     %,,,,,,,
%     Karray = round(euDist) + move;      %,,,,,,,
    Karray=[];                        %%%%%%%%

    for m = 1:N
%         for t = 1:Karray(m)  %,,,,,,,
        for t = 1:K          %%%%%%%%
            blockSize = round(MEMORY_TOTAL / N);
            counter = 0;
            while(1)
                I = abs(X' * res(:, m));
                I(m,:) = 0; 
                [~,J] = max(I, [], 1);
                S(m,nonzero(m)+1) = J;
                nonzero(m) = nonzero(m) + 1;
                counter = counter + blockSize;
                if counter >= N
                    break;
                end
            end
%             if t ~= Karray(m)   %,,,,,,,
            if t ~= K           %%%%%%%%
                    B = Xn(:, S(m, 1:nonzero(m))); 
                    res(:, m) = Xn(:, m) - B * (B \ Xn(:, m)); 
                if sum( res(:, m).^2 ) < thr
                    break;
                end 
            end
        end
    end
    for iN = 1:N	
        Val(iN,1:nonzero(iN)) = (X(:, S(iN, 1:nonzero(iN))) \ X(:, iN))';
    end   
    C = zeros(N,N);
    Ind = zeros(N,N);
    for i = 1:N
        Ind(i,1:nonzero(i)) = i;
    end
    for i = 1:N
        for j = 1:nonzero(i,:)
            C(S(i,j),Ind(i,j)) = Val(i,j);
        end
    end
end
