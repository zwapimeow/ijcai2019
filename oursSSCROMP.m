function [Val_0, Karray] = oursSSCROMP(X, K,thr)
    MEMORY_TOTAL = 0.1 * 10^9;
    [~, N] = size(X); 
    threshold = 1;
    K1 = 1;
    Xn = X; 
    Xn = cnormalize(X);
    S = zeros(N, 2*K+1); 
    Val = zeros(N, 2*K+1); 
    res = Xn; 
    nonzero = zeros(N,1);
    r = 1;
    
    euDist = X'*X;                          %,,,,,,,
    euDist = euDist + euDist';              %,,,,,,,
    euDist = sort(euDist,2,'descend');      %,,,,,,,
    euDist = euDist(:,2:K);                 %,,,,,,,
    euDist = mean(euDist,2);                %,,,,,,,
    maxDist = max(max(euDist));             %,,,,,,,
    minDist = min(euDist(euDist~=0));       %,,,,,,,
    dist = maxDist - minDist;               %,,,,,,,
    euDist = K*(euDist - minDist)/dist;     %,,,,,,,
    move = K - round(mean(euDist));         %,,,,,,,
    Karray = round(euDist) + move;          %,,,,,,,
%     Karray=[];                                %%%%%%%%

    for m = 1:N
        if (m~=1)
            [update_m,~] = find(S == m);
            if (isempty(update_m))
                X1 = X;
                Xn_1 = Xn;
            elseif (size(update_m,1)==1)
                X1 = X;
                Xn_1 = Xn;
                X1(:,update_m) =  r*((X1(:,update_m) - dot(X1(:,update_m),X1(:,m)) .* X1(:,m))/dot(X1(:,m),X1(:,m)));
                Xn_1(:,update_m) =  r*((Xn_1(:,update_m) - dot(Xn_1(:,update_m),Xn_1(:,m)) .* Xn_1(:,m))/dot(Xn_1(:,m),Xn_1(:,m)));
            elseif (size(update_m,1)>1)
                X1 = X;
                Xn_1 = Xn;
                for i = 1:size(update_m,1)
                    X1(:,update_m(i)) =  r*((X1(:,update_m(i)) - dot(X1(:,update_m(i)),X1(:,m)) .* X1(:,m))/dot(X1(:,m),X1(:,m)));
                    Xn_1(:,update_m(i)) =  r*((Xn_1(:,update_m(i)) - dot(Xn_1(:,update_m(i)),Xn_1(:,m)) .* Xn_1(:,m))/dot(Xn_1(:,m),Xn_1(:,m)));
                end   
            end
        else
            X1 = X;
            Xn_1 = Xn;
        end
        for t = 1:Karray(m)           %,,,,,,,
%         for t = 1:K                     %%%%%%%%
            blockSize = round(MEMORY_TOTAL / N);
            counter = 0;
            while(1)
                I = abs(X1' * res(:, m));
                I(m,:) = 0; 
                [~,a1_00] = max(I, [], 1);
                S(m,nonzero(m)+1) = a1_00;
                nonzero(m) = nonzero(m) + 1;
                counter = counter + blockSize;
                if counter >= N
                    break;
                end
            end
            if t ~= Karray(m)       %,,,,,,,
%             if t ~= K                 %%%%%%%%
                    B = Xn_1(:, S(m, 1:nonzero(m))); 
                    res(:, m) = Xn_1(:, m) - B * (B \ Xn_1(:, m));
                if sum( res(:, m).^2 ) < thr
                    break;
                end 
            end
        end
    end
    for iN = 1:N
        Val(iN,1:nonzero(iN)) = (X(:, S(iN, 1:nonzero(iN))) \ X(:, iN))';
    end 
    Val_0 = zeros(N,N);
    Ind = zeros(N,N);
    for i = 1:N
        Ind(i,1:nonzero(i)) = i;
    end
    for i = 1:N
        for j = 1:nonzero(i,:)
            Val_0(S(i,j),Ind(i,j)) = Val(i,j);
        end
    end
end
