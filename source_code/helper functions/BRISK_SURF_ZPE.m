% This function uses a combination of BRISK and SURF feature identification
% to create a more robust feature-based rigid registration method. Adapted
% from a Mathworks tutorial page.

%% Uses the SUP layer as the registration base

function [tformTotal, OutputView] = BRISK_SURF_ZPE(fixed_all, moving_all)

    MatchThreshold = 90;    % Works well. Increasing this increases the percent distance from a 'perfect match', i.e. more matches
    MaxRatio = 0.9;         % Works well. Increasing this increases the number of ambiguous matches, i.e. more matches
                            % 90 and 0.9 works best for OCTA
    
    ptsOriginalBRISK = detectBRISKFeatures(fixed_all);
    ptsDistortedBRISK = detectBRISKFeatures(moving_all);

    ptsOriginalSURF = detectSURFFeatures(fixed_all);
    ptsDistortedSURF = detectSURFFeatures(moving_all);

    [featuresOriginalFREAK, validPtsOriginalBRISK] = extractFeatures(fixed_all, ptsOriginalBRISK);
    [featuresDistortedFREAK, validPtsDistortedBRISK] = extractFeatures(moving_all, ptsDistortedBRISK);

    [featuresOriginalSURF, validPtsOriginalSURF] = extractFeatures(fixed_all, ptsOriginalSURF);
    [featuresDistortedSURF, validPtsDistortedSURF] = extractFeatures(moving_all, ptsDistortedSURF);

    indexPairsBRISK = matchFeatures(featuresOriginalFREAK, featuresDistortedFREAK, 'MatchThreshold', MatchThreshold, 'MaxRatio', MaxRatio, 'Unique', true);
    indexPairsSURF = matchFeatures(featuresOriginalSURF, featuresDistortedSURF);

    matchedOriginalBRISK = validPtsOriginalBRISK(indexPairsBRISK(:,1));
    matchedDistortedBRISK = validPtsDistortedBRISK(indexPairsBRISK(:,2));

    matchedOriginalSURF = validPtsOriginalSURF(indexPairsSURF(:,1));
    matchedDistortedSURF = validPtsDistortedSURF(indexPairsSURF(:,2));

    matchedOriginalXY = [matchedOriginalSURF.Location; matchedOriginalBRISK.Location];
    matchedDistortedXY = [matchedDistortedSURF.Location; matchedDistortedBRISK.Location];
      
    % If there aren't enough matches, discard the strip
    if (size(matchedOriginalXY,1) < 3) || (size(matchedDistortedXY,1) < 3)
        disp('SURF Registration failed: Not enough keypoint matches')
        return
    end

    [tformTotal, inlierDistortedXY, inlierOriginalXY] = estimateGeometricTransform(matchedDistortedXY, matchedOriginalXY, 'affine', 'MaxNumTrials', 50000);
    OutputView = imref2d(size(fixed_all));

    %% If the shear of the transformation is greater than 0.3 in either X or Y, registration failed
    if abs(tformTotal.T(1,2)) > 0.3 || abs(tformTotal.T(2,1)) > 0.3
        disp('BRISK/SURF Registration failed: Sheer transformation too high')
        return
    end
end