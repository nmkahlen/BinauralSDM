% Copyright (c) Facebook, Inc. and its affiliates.

function [BRIR_DSER, BRIR_LR, BRIR_DS, BRIR_ER] = Normalize(BRIR_data, BRIR_DSER, BRIR_LR, BRIR_DS, BRIR_ER)
% This function normalizes the response to the frontal direction,
% either to direct sound, direct sound and early reflections, or overall. Always done to frontal
% directions. A-weighting is used.

%   Author: Nils Meyer-Kahlen
%   Last modified: 11/09/2022


% use response closest to frontal

if ~strcmp(BRIR_data.NormalizationMode, 'none')

    [x, y, z] = sph2cart(BRIR_data.Directions(:, 1) * pi / 180 , ...
        BRIR_data.Directions(:, 2) * pi / 180, 1)
    idxFront  = knnsearch( [x,y,z], [1, 0, 0]);

    switch BRIR_data.NormalizationMode
        case 'full'

            BRIR_selected = BRIR_LR;
            BRIR_selected(1:size(BRIR_DSER, 1, 1), :) = ...
                BRIR_selected(1:size(BRIR_DSER, 1, 1), :) + BRIR_DSER(:, :, idxZero);

        case 'direct'

            % normalize to direct level,
            BRIR_selected = BRIR_DS(:, :, idxFront);

        case 'direct+early'

            % normalize to direct and ES sound level,
            % just to match Hannes' test
            BRIR_selected = BRIR_DSER(:, :, idxFront);

    end

    % A weighting
    A_filter = weightingFilter('A-weighting',  BRIR_data.fs);
    BRIR_selected_A_weighted = A_filter(BRIR_selected);

    % mean RMS
    levelA = mean(rms(BRIR_selected_A_weighted));

    BRIR_DS = BRIR_DS / levelA;
    BRIR_ER = BRIR_ER / levelA;
    BRIR_LR = BRIR_LR / levelA;
    BRIR_DSER = BRIR_DSER / levelA;

    disp(['Normalized to ' BRIR_data.NormalizationMode '. Changed by ' num2str(db(levelA)), ' dB'])

end

end




