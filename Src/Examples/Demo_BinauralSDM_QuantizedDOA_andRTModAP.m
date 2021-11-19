% Copyright (c) Facebook, Inc. and its affiliates.

%% Description

% This script is an example showing how to use the functions provided in
% the BinauralSDM repository. The BRIRs generated by this script include
% the steps described in Amengual et al. 2020 - DOA Postprocessing 
% (smoothing and quantization) and RTMod+AP equalization.

% The overall process is as follows:
% - Spatial data from a multichannel RIR is obtained using the SDM
% utilizing the functions from the SDM Toolbox (Tervo & Patynen).
% - The spatial information is smoothed and quantized as proposed in
% Amengual et al. 2020.
% - Preliminary BRIRs are synthesized by selecting the closest directions
% from an HRIR dataset. In this example, a dataset of KU100 dummy head is
% utilized. The dataset was generated by Bernschutz et al. from the Audio
% Group of TH Cologne. The HRIR dataset is downloaded on the fly during the
% example.
% - The reverberation time of the preliminary BRIRs is corrected using the
% RTMod method introduced in Amengual et al. 2020.
% - After RTMod correction, a cascade of 3 all-pass filters is applied to
% both left and right channels to improve the diffuseness of the late
% reverberation tail. This process is also introduced in Amengual et al.
% 2020.
% - Responses for arbitrarily defined head orientations are generated
% using the previous steps by rotating the DOA vectors.
% - The data is saved in a user defined folder.
%
%
% References: 
% - (Tervo et al. 2013) - "Spatial Decomposition Method for Room Impulse
% Responses", JAES, 2013.
% - (SDM Toolbox) - "SDM Toolbox", S. Tervo and J. Patynen, Mathworks
% Exchange
% - (Amengual et al. 2020) - "Optimizations of the Spatial Decomposition 
% Method for Binaural Reproduction", JAES 2020.

% Author: Sebastia Amengual (samengual@fb.com)
% Last modified: 11/17/2021

clear; clc; close all;

% set current folder to the one containing this script for relative paths 
% to work properly
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clear tmp;

% add base folder to path for dependencies to be available
addpath(genpath('../../'));

%% Analysis parameter initialization
time_start = tic; % start measuring execution time

% Analysis parameters
MicArray        = 'FRL_10cm';           % FRL Array is 7 mics, 10cm diameter, with central sensor. Supported geometries are FRL_10cm, FRL_5cm, 
                                        % Tetramic and Eigenmike. Modify the file create_MicGeometry to add other geometries (or contact us, we'd be happy to help).
Room            = 'ExampleRoom';        % Name of the room. RIR file name must follow the convention RoomName_SX_RX.wav
SourcePos       = 'S1';                 % Source Position. RIR file name must follow the convention RoomName_SX_RX.wav
ReceiverPos     = 'R1';                 % Receiver Position. RIR file name must follow the convention RoomName_SX_RX.wav
Database_Path   = '../../Data/RIRs/';   % Relative path to folder containing the multichannel RIR
fs              = 48e3;                 % Sampling Rate (in Hz). Only 48 kHz is recommended. Other sampling rates have not been tested.
MixingTime      = 0.08;                 % Mixing time (in seconds) of the room for rendering. Data after the mixing time will be rendered 
                                        % as a single direction independent reverb tail and AP rendering will be applied.
DOASmooth       = 16;                   % Window length (in samples) for smoothing of DOA information. 16 samples is a good compromise for noise 
                                        % reduction and time resolution.
DOAOnsetLength  = 128;                  % Length (in samples) for assignment of a constant (averaged) DOA for the onset / direct sound.
DenoiseFlag     = true;                 % Flag to perform noise floor compensation on the multichannel RIR. If set, it ensures that the RIR decays 
                                        % progressively and removes rendering artifacts due to high noise floor in the RIR.
FilterRawFlag   = true;                 % Flag to perform band pass filtering on the multichannel RIR prior to DOA estimation. If set, only
                                        % information between 200 Hz and 8 kHz (by default) will be used for DOA estimation. This helps increasing 
                                        % robustness of the estimation. See create_BRIR_data.m for customization of the filtering.
AlignDOAFlag    = true;                 % If this flag is set, the DOA data will be rotated so the direct sound is aligned to 0,0 (az, el).
SpeedSound      = 345;                  % Speed of sound in m/s (for SDM Toolbox DOA analysis)
WinLen          = 62;                   % Window Length (in samples) for SDM DOA analysis. For fs = 48kHz, sizes between 36 and 64 seem appropriate. 
                                        % The optimal size might be room dependent. See Tervo et al. 2013 and Amengual et al. 2020 for a discussion.

% Initialize SRIR data struct
SRIR_data = create_SRIR_data('MicArray',MicArray,...
                             'Room',Room,...
                             'SourcePos',SourcePos,...
                             'ReceiverPos',ReceiverPos,...
                             'Database_Path',Database_Path,...
                             'fs',fs,...
                             'MixingTime',MixingTime,...
                             'DOASmooth',DOASmooth,...
                             'DOAOnsetLength',DOAOnsetLength,...
                             'Denoise',DenoiseFlag,...
                             'FilterRaw',FilterRawFlag,...
                             'AlignDOA',AlignDOAFlag);
clear MicArray Room SourcePos ReceiverPos Database_Path fs MixingTime ...
    DOASmooth BRIRLength DenoiseFlag FilterRawFlag AlignDOAFlag;

% Initialize SDM analysis struct (from SDM Toolbox)
SDM_Struct = createSDMStruct('c',SpeedSound,...
                             'fs',SRIR_data.fs,...
                             'micLocs',SRIR_data.ArrayGeometry,...
                             'winLen',WinLen);
clear SpeedSound WinLen;

%% Download a HRIR from the Cologne audio team server
% (skipped automatically if the HRIR dataset already exists)
HRIR_URL = 'https://zenodo.org/record/3928297/files/HRIR_FULL2DEG.sofa?download=1';
HRIR_Folder   = '../../Data/HRIRs/';
HRIR_Filename = 'KU100_HRIR_FULL2DEG_Koeln.sofa';
HRIR_Subject  = 'KU100';    % Name of the HRIR subject (only used for naming purposes while saving).
HRIR_Type     = 'SOFA';     % File format of the HRIR. Only SOFA is supported for now.

[~,~] = mkdir(HRIR_Folder); % ignore warning if directory already exists
HRIR_Path = fullfile(HRIR_Folder, HRIR_Filename);
clear HRIR_Folder HRIR_Filename;

fprintf('Downloading HRIR dataset ... ');
if isfile(HRIR_Path)
    fprintf('skipped.\n\n');
else
    websave(HRIR_Path, HRIR_URL);
    fprintf('done.\n\n');
end; clear HRIR_URL;

%% Rendering parameters
DOAAzOffset      = 0.0;             % Azimuth rotation in degrees aplied after DOA estmation and before DOA quantization / BRIR rendering.
DOAElOffset      = 0.0;             % Elevation rotation in degrees aplied after DOA estmation and before DOA quantization / BRIR rendering. 
QuantizeDOAFlag = true;             % Flag to determine if DOA information must me quantized.
DOADirections   = 50;               % Number of directions to which the spatial information will be quantized using a Lebedev grid.
BandsPerOctave  = 3;                % Bands per octave used for RT60 equalization
EqTxx           = 30;               % Txx used for RT60 equalization. For very small rooms or low SNR measurements, T20 is recommended. Otherwise, T30 is recommended.
RTModRegFreq    = false;            % Regularization frequncy in Hz above which the RTmod modification of the late reverberation should be regularized.
NamingCond      = sprintf('Quantized%dDOA', DOADirections); % String used for naming purposes, useful when rendering variations of the sasme RIR.
BRIR_Length     = 0.0;              % Length of BRIR in seconds, will be chosen by analysis of the room reverberation time if unspecified.
BRIRAtten       = 30;               % Attenuation of the rendered BRIRs (in dB). Useful to avoid clipping. Use the same value when rendering various
                                    % positions in the same room to maintain relative level differences.
AzOrient        = (-180:2:180)';    % Render BRIRs every 2 degrees in azimuth.
ElOrient        = (-90:5:90)';      % Render BRIRs every 5 degrees in elevation.
PlotAnalysisFlag = true;            % Flag to determine if analysis results should be plotted.ahaha
PlotExportFlag   = true;            % Flag to determine if analysis result plots should be exported.
DestinationPath = '../../Data/RenderedBRIRs/'; % Folder where the resulting BRIRs will be saved.

% Append configuration to destination path
DestinationPath = fullfile(DestinationPath, ...
    strrep(HRIR_Subject, ' ', '_'), ...
    sprintf('%s_%s_%s', SRIR_data.Room, SRIR_data.SourcePos, SRIR_data.ReceiverPos), ...
    NamingCond);

% Initialize re-synthesized BRIR struct
BRIR_data = create_BRIR_data('MixingTime',SRIR_data.MixingTime,...
                             'HRTF_Subject',HRIR_Subject,...
                             'HRTF_Type',HRIR_Type,...
                             'HRTF_Path',HRIR_Path,...
                             'BandsPerOctave',BandsPerOctave,...
                             'EqTxx',EqTxx,...
                             'RTModRegFreq',RTModRegFreq,...
                             'Length',BRIR_Length,...
                             'RenderingCondition',NamingCond,...
                             'Attenuation',BRIRAtten,...
                             'AzOrient',AzOrient,...
                             'ElOrient',ElOrient,...
                             'DOAAzOffset',DOAAzOffset,...
                             'DOAElOffset',DOAElOffset,...
                             'QuantizeDOAFlag',QuantizeDOAFlag,...
                             'DOADirections',DOADirections,...
                             'DestinationPath',DestinationPath,...
                             'fs',SRIR_data.fs);
clear HRIR_Subject HRIR_Type HRIR_Path BandsPerOctave EqTxx RTModRegFreq ...
    NamingCond BRIRAtten AzOrient ElOrient ...
    DOAAzOffset DOAElOffset QuantizeDOAFlag DOADirections DestinationPath;

% Initialize visualization struct (from SDM Toolbox, with additions)
Plot_data = createVisualizationStruct(...
    'fs', SRIR_data.fs, 'DefaultRoom', 'Small', ... % or 'VerySmall, 'Medium', 'Large'
    'name', strrep(sprintf('%s_%s_%s', ...
    SRIR_data.Room, SRIR_data.SourcePos, SRIR_data.ReceiverPos), '_', '\_'));
Plot_data.PlotAnalysisFlag = PlotAnalysisFlag;
Plot_data.PlotExportFlag = PlotExportFlag;
Plot_data.DestinationPath = BRIR_data.DestinationPath;
clear PlotAnalysisFlag PlotExportFlag;

%% Analysis
% Estimate directional information using SDM. This function is a wrapper of
% the SDM Toolbox DOA estimation (using TDOA analysis) to include some 
% post-processing. The actual DOA estimation is performed by the SDMPar.m 
% function of the SDM Toolbox.
SRIR_data = Analyze_SRIR(SRIR_data, SDM_Struct);

if Plot_data.PlotAnalysisFlag
    % Plot analysis results
    plot_name = [Plot_data.name, '_spatio_temporal'];
    if SRIR_data.AlignDOA; plot_name = [plot_name, '_aligned']; end
    
    Plot_Spec(SRIR_data, Plot_data, [Plot_data.name, '_time_frequency']);
    Plot_DOA(SRIR_data, Plot_data, plot_name);
    clear plot_name;
end

%% Synthesis
% 1. Pre-processing operations (massage HRTF directions, resolve DOA NaNs).

% Read HRTF dataset for re-synthesis
HRTF_data = Read_HRTF(BRIR_data);

[SRIR_data, BRIR_data, HRTF_data, HRTF_TransL, HRTF_TransR] = ...
    PreProcess_Synthesize_SDM_Binaural(SRIR_data, BRIR_data, HRTF_data);

% -----------------------------------------------------------------------
% 2. Rotate DOA information, if required

if BRIR_data.DOAAzOffset || BRIR_data.DOAElOffset
    SRIR_data = Rotate_DOA(SRIR_data, ...
        BRIR_data.DOAAzOffset, BRIR_data.DOAElOffset);

    if Plot_data.PlotAnalysisFlag
        Plot_DOA(SRIR_data, Plot_data, [Plot_data.name, '_spatio_temporal_rotated']);
    end
end

% -----------------------------------------------------------------------
% 3. Quantize DOA information, if required

if BRIR_data.QuantizeDOAFlag
    [SRIR_data, ~] = QuantizeDOA(SRIR_data, ...
        BRIR_data.DOADirections, SRIR_data.DOAOnsetLength);
    
    if Plot_data.PlotAnalysisFlag
        Plot_DOA(SRIR_data, Plot_data, [Plot_data.name, '_spatio_temporal_quantized']);
    end
end

% -----------------------------------------------------------------------
% 4. Compute parameters for RTMod Compensation

% Synthesize one direction to extract the reverb compensation - solving the
% SDM synthesis spectral whitening
BRIR_Pre = Synthesize_SDM_Binaural(...
    SRIR_data, BRIR_data, HRTF_TransL, HRTF_TransR, [0, 0], true);

% Using the pressure RIR as a reference for the reverberation compensation
BRIR_data.ReferenceBRIR = [SRIR_data.P_RIR, SRIR_data.P_RIR];

% Get the desired T30 from the Pressure RIR and the actual T30 from one
% rendered BRIR
[BRIR_data.DesiredT30, BRIR_data.OriginalT30, BRIR_data.RTFreqVector] = ...
    GetReverbTime(SRIR_data, BRIR_Pre, BRIR_data.BandsPerOctave, BRIR_data.EqTxx);
clear BRIR_Pre;

% -----------------------------------------------------------------------
% 5. Render BRIRs with RTMod compensation for the specified directions

nDirs = length(BRIR_data.Directions);

% Render early reflections
hbar = parfor_progressbar(nDirs, 'Please wait, rendering (step 1/2) ...');
parfor iDir = 1 : nDirs
    hbar.iterate(); %#ok<PFBNS>
    BRIR_TimeDataTemp = Synthesize_SDM_Binaural(SRIR_data, BRIR_data, ...
        HRTF_TransL, HRTF_TransR, BRIR_data.Directions(iDir, :), false);
    BRIR_TimeDataTemp = Modify_Reverb_Slope(BRIR_data, BRIR_TimeDataTemp);
    BRIR_TimeData(:, :, iDir) = BRIR_TimeDataTemp;
end
close(hbar);
clear iDir hbar;

% Render late reverb
BRIR_full = Synthesize_SDM_Binaural(SRIR_data, BRIR_data, ...
    HRTF_TransL, HRTF_TransR, [0, 0], true);
BRIR_full = Modify_Reverb_Slope(BRIR_data, BRIR_full, Plot_data);

% Remove leading zeros
[BRIR_TimeData, BRIR_full] = Remove_BRIR_Delay(BRIR_TimeData, BRIR_full, -20);

% Split the BRIR
[BRIR_DSER, BRIR_LR, BRIR_DS, BRIR_ER] = Split_BRIR(...
    BRIR_TimeData, BRIR_full, BRIR_data.MixingTime, BRIR_data.fs, 256);
clear BRIR_full;

% -----------------------------------------------------------------------
% 6. Apply AP processing for the late reverb

% AllPass filtering for the late reverb (increasing diffuseness and
% smoothing out the EDC)
BRIR_data.allpass_delays = [37, 113, 215, 347]; % in samples
BRIR_data.allpass_RT = [0.1, 0.1, 0.1, 0.1];    % in seconds

for iAllPass = 1:3
    BRIR_LR(:,1) = allpass_filter(BRIR_LR(:,1), ...
        BRIR_data.allpass_delays(iAllPass), 0.1, BRIR_data.fs);
    BRIR_LR(:,2) = allpass_filter(BRIR_LR(:,2), ...
        BRIR_data.allpass_delays(iAllPass), 0.1, BRIR_data.fs);
end; clear iAllPass;

% -----------------------------------------------------------------------
% 7. Save BRIRs

if Plot_data.PlotAnalysisFlag
    Plot_BRIR(BRIR_data, BRIR_DS, BRIR_ER, BRIR_LR, Plot_data);
end

hbar = parfor_progressbar(nDirs + 1, 'Please wait, saving (step 2/2) ...');
for iDir = 1 : nDirs
    hbar.iterate();
    if iDir == 1 % export identical late reverberation only once
        BRIR_LR_export = BRIR_LR;
    else
        BRIR_LR_export = [];
    end
    SaveBRIR(BRIR_data, BRIR_DS(:, :, iDir), BRIR_DSER(:, :, iDir), ...
        BRIR_ER(:, :, iDir), BRIR_LR_export, BRIR_data.Directions(iDir, :));
end
hbar.iterate();
SaveRenderingStructs(SRIR_data, BRIR_data);
close(hbar);
clear iDir hbar BRIR_LR_export nDirs;

%%
time_exec = toc(time_start);
fprintf('\n... completed in %.0fh %.0fm %.0fs.\n', ...
    time_exec/60/60, time_exec/60, mod(time_exec, 60));
