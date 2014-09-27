function [segQuality connList avgLik avgSP]= evalSeg(im_test, im_gt, seg_fgCurrent, seg_fgPrior)

load topologyPrior;

% Handle the case of Empty Segments ...
if( numel( unique(seg_fgCurrent) ) == 1 || numel( unique(seg_fgPrior) ) == 1 )
    segQuality = -log(eps);
    connList = [];
    avgLik = 0;
    avgSP = 0;
    return;
end

% Initialize Parameters
    
    avgLik = 0;
    avgSP = 0;
    spPrior = 0;

% Connected Component Analysis
    CCPrior   = bwconncomp(seg_fgPrior);
    CCCurrent = bwconncomp(seg_fgCurrent);
    for priorIter = 1:CCPrior.NumObjects
        CCPrior.numPixels{priorIter} = numel(CCPrior.PixelIdxList{priorIter});
    end
    for currIter = 1:CCCurrent.NumObjects
        CCCurrent.numPixels{currIter} = numel(CCCurrent.PixelIdxList{currIter});
    end
    connList = cell(CCPrior.NumObjects, 1);
    connStat = cell(CCPrior.NumObjects, 1);
    currObjs = ones(CCCurrent.NumObjects, 1);
    dataLik = cell(CCPrior.NumObjects,1);
% Actual Analysis Routine
    for priorIter = 1:CCPrior.NumObjects
        I1 = zeros(size(seg_fgPrior));
        I2 = zeros(size(seg_fgCurrent));
        priorContour = CCPrior.PixelIdxList;
        I1(priorContour{priorIter}(:)) = 1;
        for currIter = 1:CCCurrent.NumObjects
            Itemp = zeros(size(seg_fgPrior));
            currContour = CCCurrent.PixelIdxList;
            Itemp(currContour{currIter}(:)) = 1;
            if( sum( im2bw(I1(:)) & im2bw(Itemp(:)) ) > 0 )
                connList{priorIter} = [connList{priorIter}; currIter];
                connStat{priorIter} = [connStat{priorIter};
                    CCCurrent.numPixels{currIter}/ CCPrior.numPixels{priorIter}];
                I2(currContour{currIter}(:)) = 1;
                currObjs(currIter) = 0;
            end
        end
        if(sum(I2(:)) == 0 ) % Modeling object disappearance
            connList{priorIter} = [];
            connStat{priorIter} = 0;
        end
        [fg_proj priorHist] = vrl_msPatchHist(im_gt , find(I1~=0), 16);
        [fg_proj currHist] = vrl_msPatchHist(im_test , find(I2~=0), 16);
        dataLik{priorIter} = [dataLik{priorIter}; sum( min( priorHist(:), currHist(:) ) )];
        % Detect a split for sure
        if( numel(connList{priorIter}) > 1)
            spPrior = normpdf( sum(connStat{priorIter})*100, learntParam.splitSum.mean, learntParam.splitSum.std/3);         
            display('SPLIT');            
        elseif( sum(I2(:) ~= 0 ) < sum(I1(:)~=0) )
            % This denotes a shrinkage
            connStat{priorIter}*100
            spPrior = normpdf(connStat{priorIter}*100, learntParam.shrinkage.mean, learntParam.shrinkage.std/3);
            display('SHRINK');
        elseif( sum(I2(:) ~= 0 ) > sum(I1(:)~=0) )
            if( sum(I2(:) ~= 0) < 1.25 * sum(I1(:)~=0)  )
                connStat{priorIter}*100
                spPrior = normpdf(connStat{priorIter}*100, learntParam.shrinkage.mean+20, learntParam.shrinkage.std/3);
                display('EXPAND');
            else
                spPrior = normpdf(connStat{priorIter}*100, learntParam.merge.mean+50, learntParam.merge.std/3);
                display('MERGE');
            end
        end
        segIntermediate(priorIter) = spPrior * dataLik{priorIter};
        avgLik = avgLik + dataLik{priorIter};
        avgSP = avgSP + spPrior;
    end
    avgLik = avgLik / CCPrior.NumObjects;
    avgSP = avgSP / CCPrior.NumObjects;
    segQuality = -log( mean(segIntermediate) * double(sum(currObjs) == 0) + eps );