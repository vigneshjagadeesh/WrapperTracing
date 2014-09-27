%% Inference specific Computations
if( inference == 1)
    [minQual currTheta] = min(qualityMat(:));
    thetaOpt = [thetaOpt; theta( currTheta, : )];
    segFgPrior = currProposals{currTheta};
elseif( inference == 2)
    %% Iterate over all states in this time instant
    for stateIter = 1:numParams
        [max_val psiDP(stateIter, frameIter) ] = max( delDP(:, frameIter-1) .* qualityMat(stateIter,:)' ); % psi simply stores the indices ...
        delDP(stateIter, frameIter) = obsDP(stateIter, frameIter) .* max_val;
    end
    [optimal_value currTheta ] = max( delDP(:,frameIter) );
    thetaOpt = [thetaOpt; theta( currTheta, : )];
    segFgPrior = currProposals{currTheta};
elseif( inference == 3)
    [minQual currTheta] = min(qualityMat(:));
    thetaOpt = [thetaOpt; theta( currTheta, : )];
    segFgPrior = currProposals{currTheta};
    wPF = wPF .* exp( - qualityMat );
    wPF = wPF ./ sum(wPF);
    % Resample Particles
    resamp = randsample(1:numParams, numParams, true, wPF);
    theta = theta( resamp, :);
    theta = randn( numParams, 2) * diag(tvPF) + theta;
end
currOptima(:,:,frameIter) = segFgPrior; % Accumulate Results for 3D Recon