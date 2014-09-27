function [fgPost bgPost fgndFrac] = refineContourPrior( imCurr, imPrev, segFgPrior, max_dist, USEFLOW )


%% Compute Likelihood Potential using Multi Scale Histograms
Ifg = find(segFgPrior ~= 0);
Ibg =  find(bwdist(segFgPrior) < sqrt(numel(Ifg)/pi) & bwdist(segFgPrior) > 0); % find(seg_fgPrior == 0); %
[fg_proj fg_hist] = vrl_msPatchHist(imPrev, Ifg, 8);
fgndProject = vrl_msPatchHistbp(imCurr, fg_hist);
[bg_proj bg_hist] = vrl_msPatchHist(imPrev, Ibg, 8);
bgndProject = vrl_msPatchHistbp(imCurr, bg_hist);

[Vx,Vy,reliab] = optFlowLk( imPrev, imCurr, [], 2, 1.2, 3e-6, 0 );
connComps = bwconncomp(segFgPrior);
for connCompIter = 1:connComps.NumObjects
    currCC = zeros(size(imCurr));
    currCC(connComps.PixelIdxList{connCompIter}) = 1;
    [signed_dist contour_band] = compute_signed_dist(currCC, 10);
    fgndFrac(connCompIter) = sum(currCC==1) / (sum(currCC==1) + sum(currCC==0) );
    [S Sx Sy] = f_GennerateSigDist(currCC);
    [F Fx Fy] = f_GennerateSigDist( Sx .* Vx + Sy .* Vy );
    S = S - F * USEFLOW ; %% Contour Deformation Induced by the Motion Prior
    currContour = S > 0;
    expandVar = ( (1-fgndFrac(connCompIter)) * max_dist ) / 30; % Try Learning this Offline
    fgPrior(:,:,connCompIter) = exp( - bwdist(currContour) / expandVar ) ;
    bgP = 1 - fgPrior(:,:,connCompIter) ;
    bgP(bgP == 0) = min(bgP(bgP~=0));
    bgPrior(:,:,connCompIter) = bgP;
    % figure(100); imshow(im_test,[]);  hold on; contour(seg_fgPrior, [0 0], 'color', 'green'); contour(S, [0 0], 'color', 'red');  quiver(Vx, Vy); hold off;
end


%% Final Likelihoods
fgPrior = adjustRange( max(fgPrior, [], 3), [0 1], [0 1]);
bgPrior = adjustRange( min(bgPrior, [], 3), [0 1], [0 1]);
fgPost = -log(fgndProject .*fgPrior  + 0.05);
bgPost = -log(bgndProject .*bgPrior  + 0.05);
segML = fgPost>bgPost;