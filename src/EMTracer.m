% Initialize the Tracing Procedure ... set all  parameters in the script
% below
initEMTracer;

for frameIter = 2:numSlices-1
    %% Change variables over time
       imPrev = I_stack(:,:,frameIter);
       imCurr = I_stack(:,:,frameIter+1);     
       
   %% SRS-DP
    [fgPost bgPost] = contourSearch( imCurr, imPrev, segFgPrior, max_dist, USEFLOW );
       
    %% Interaction Potentials ...
       imSmooth = imfilter(imCurr, fspecial('gauss', [7 7], 1.2));
       % imSmooth = kovesi_anisodiff(double(imCurr), 40, 20, .25, 1);
       im_test_smooth = double(reshape(imSmooth, noRows*noCols, 1));       
       % visOutputs(fgndProject, bgndProject, fgPrior, bgPrior, fgPost, bgPost, imCurr, imSmooth, imPrev, segML, segFgPrior, Vx, Vy, nlinkRange(1), [1 1 0 1]);  return;       
       for thetaIter = 1:numParams
           [weights,w_dist] = vrl_edgeweight(edges, im_test_smooth', [noRows noCols], theta(thetaIter, 1) ); 
           % weights = vrlConstrastPotentials(imCurr, edges, theta(thetaIter, 1) );
           currProposals{thetaIter} =  reshape( vrl_gc( int32([2, noRows * noCols, no_nlinks]), ...                     
                                                        double([fgPost(:) bgPost(:)]'), ...        
                                                        [edges-1, theta(thetaIter, 2) * weights]' ...       
                                                        ) , noRows, noCols );   
           plotNum = mod(thetaIter, 9); if(plotNum == 0), plotNum = 9; end
           figure( ceil(thetaIter/9) ); subplot(3,3, plotNum); imshow(imCurr); hold on; contour( currProposals{thetaIter}-0.5, [0 0], 'r'); hold off;
           if( inference == 1 || inference == 3)
                   qualityMat(thetaIter) = scoreTransitions(imCurr, imPrev, currProposals{thetaIter}, segFgPrior);
           elseif( inference == 2 )  % segHMMPrior is the Transition Matrix
                   for prevIter = 1:numParams
                        qualityMat(thetaIter, prevIter) = exp( - evalSeg(imCurr, imPrev, currProposals{thetaIter}, prevProposals{thetaIter}) );
                   end
           end
       end       
       regionSearch;

    if(numel(unique(segFgPrior)) == 1 || numel(qualityMat) == 1), display('END OF STACK REACHED --- QUITTING'); break; end
    figure(3); imshow(imCurr); hold on; contour(currOptima(:,:,frameIter), [0 0], 'g'); hold off;
end