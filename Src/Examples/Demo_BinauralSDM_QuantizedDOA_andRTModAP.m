% Copyright (c) Facebook, Inc. and its affiliates.

%% Description

% This script is an example showing how to use the functions provided in
% the BinauralSDM repository. The BRIRs generated by this script include
% the steps described in Amengual et al. 2020 - DOA Postprocessing 
% (smoothing and quantization) and RTMod+AP equalization. 
% 
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
% References:
% - (Tervo et al. 2013) - "Spatial Decomposition Method for Room Impulse
% Responses", JAES, 2013.
% - (SDM Toolbox) - "SDM Toolbox", S. Tervo and J. Patynen, Mathworks
% Exchange
% - (Amengual et al. 2020) - "Optimizations of the Spatial Decomposition 
% Method for Binaural Reproduction", JAES 2020.
% 
% Author: Sebastia Amengual (samengual@fb.com)
% Last modified: 11/17/2021
% 
% IMPORTANT NOTE: To ensure that relative paths work, navigate to the
% folder containing this script before executing.

%% Analysis parameter initialization

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
BRIRLength      = 0.5;                  % Duration of the rendered BRIRs (in seconds)
DenoiseFlag     = 1;                    % Flag to perform noise floor compensation on the multichannel RIR. This ensures that the RIR decays 
                                        % progressively and removes rendering artifacts due to high noise floor in the RIR.
FilterRawFlag   = 1;                    % Flag to perform band pass filtering on the multichannel RIR prior to DOA estimation. If active, only
                                        % information between 200Hz and 8kHz (by default) will be used for DOA estimation. This helps increasing 
                                        % robustness of the estimation. See create_BRIR_data.m for customization of the filtering.
AlignDOA        = 1;                    % If this flag is set to 1, the DOA data will be rotated so the direct sound is aligned to 0,0 (az, el).
SpeedSound      = 345;                  % Speed of sound in m/s (for SDM Toolbox DOA analysis)
WinLen          = 62;                   % Window Length (in samples) for SDM DOA analysis. For fs = 48kHz, sizes between 36 and 64 seem appropriate. 
                                        % The optimal size might be room dependent. See Tervo et al. 2013 and Amengual et al. 2020 for a discussion.

% Initialize SRIR data struct
SRIR_data = create_SRIR_data('MicArray', MicArray,...
                             'Room',Room,...
                             'SourcePos',SourcePos,...
                             'ReceiverPos',ReceiverPos,...
                             'Database_Path',Database_Path,...
                             'fs',fs,...
                             'MixingTime',MixingTime,...
                             'DOASmooth',DOASmooth,...
                             'Length',BRIRLength,...
                             'Denoise',DenoiseFlag,...
                             'FilterRaw',FilterRawFlag,...
                             'AlignDOA',AlignDOA);

% Initialize SDM analysis struct (from SDM Toolbox)
SDM_Struct = createSDMStruct('c',SpeedSound,...
                             'fs',SRIR_data.fs,...
                             'micLocs',SRIR_data.ArrayGeometry,...
                             'winLen',62);

%% Download a HRIR from the Cologne audio team server
% This step can be skipped if you already have a HRIR dataset

HRIRurl = 'https://zenodo.org/record/3928297/files/HRIR_FULL2DEG.sofa?download=1';
HRIRfolder = '../../Data/HRIRs/';
HRIRfilename = 'KU100_HRIR_FULL2DEG_Koeln.sofa';
disp('Downloading HRIR Dataset');
if ~exist(HRIRfolder,'dir')
    mkdir(HRIRfolder)
end
HRIRsave = websave([HRIRfolder HRIRfilename],HRIRurl);

%% Rendering parameters
QuantizeDOAFlag     = 1;                % Flag to determine if DOA information must me quantized.
DOADirections   = 50;               % Number of directions to which the spatial information will be quantized using a Lebedev grid
HRIR_Subject    = 'KU100';          % Name of the HRIR subject (only used for naming purposes while saving).
HRIR_Type       = 'SOFA';           % File format of the HRIR. Only SOFA is supported for now.
HRIR_Path       = '../../Data/HRIRs/KU100_HRIR_FULL2DEG_Koeln.sofa'; % Relative path to the HRIR.
NamingCond      = ['Quantized' num2str(DOADirections) 'DOA']; % String used for naming purposes, useful when rendering variations of the same RIR.
BRIRAtten       = 30;               % Attenuation of the rendered BRIRs (in dB). Useful to avoid clipping. Use the same value when rendering various
                                    % positions in the same room to maintain relative level differences.
AzOrient        = (-180:2:180)';    % Render BRIRs every 2 degrees in azimuth
ElOrient        = (-90:5:90)';      % Render BRIRs every 5 degrees in elevation
DestinationPath = 'C:/Data/RenderedBRIRs/'; % Folder where the resulting BRIRs will be saved
BandsPerOctave  = 1;                % Bands per octave used for RT60 equalization
EqTxx           = 20;                % Txx used for RT60 equalization. For very small rooms or low SNR measurements, T20 is recommended. Otherwise, T30 is recommended.

% Initialize re-synthesized BRIR struct
BRIR_data = create_BRIR_data('MixingTime',MixingTime,...
                             'HRTF_Subject',HRIR_Subject,...
                             'HRTF_Type',HRIR_Type,...
                             'HRTF_Path',HRIR_Path,...
                             'Length',BRIRLength,...
                             'RenderingCondition',NamingCond,...
                             'Attenuation',BRIRAtten,...
                             'AzOrient',AzOrient,...
                             'ElOrient',ElOrient,...
                             'QuantizeDOAFlag',QuantizeDOAFlag,...
                             'DOADirections',DOADirections,...
                             'DestinationPath',DestinationPath,...
                             'fs',fs,...
                             'BandsPerOctave',1,...
                             'EqTxx',20);

% Read HRTF dataset for re-synthesis
HRTF = Read_HRTF(BRIR_data);

%% Analysis

% Estimate directional information using SDM. This function is a wrapper of
% the SDM Toolbox DOA estimation (using TDOA analysis) to include some 
% post-processing. The actual DOA estimation is performed by the SDMPar.m 
% function of the SDM Toolbox.
SRIR_data = Analyze_SRIR(SRIR_data, SDM_Struct); 
                                                

%% Synthesis

% 1. Pre-processing operations (massage HRTF directions, resolve DOA NaNs).

[SRIR_data, BRIR_data, HRTF_data, HRTF_TransL, HRTF_TransR] = PreProcess_Synthesize_SDM_Binaural(SRIR_data, BRIR_data, HRTF);

% -----------------------------------------------------------------------
%%% 2. Quantize DOA information, if required

if BRIR_data.QuantizeDOAFlag == 1
    [SRIR_data, idx] = QuantizeDOA(SRIR_data, BRIR_data.DOADirections, 128);
end

% -----------------------------------------------------------------------
%%% 3. Compute parameters for RTMod Compensation

% Synthesize one direction to extract the reverb compensation - solving the
% SDM synthesis spectral whitening
PreBRIR = Synthesize_SDM_Binaural(SRIR_data, BRIR_data, HRTF_TransL, HRTF_TransR, [0 0],1);

% Using the pressure RIR as a reference for the reverberation compensation
BRIR_data.ReferenceBRIR = [SRIR_data.P_RIR SRIR_data.P_RIR];

% Get the desired T30 from the Pressure RIR and the actual T30 from one
% rendered BRIR
[DesiredT30, OriginalT30, RTFreqVector] = GetReverbTime(SRIR_data, PreBRIR,BRIR_data.BandsPerOctave,BRIR_data.EqTxx); 

% -----------------------------------------------------------------------
% 4. Render BRIRs with RTMod compensation for the specified directions

% Initialize BRIR matrix
BRIR_Early = zeros((BRIR_data.MixingTime+BRIR_data.TimeGuard)*BRIR_data.fs,2,length(BRIR_data.Directions));

% Render BRIRs
nDirs = length(BRIR_data.Directions);

% Render early reflections
hbar = parfor_progressbar(nDirs,'Please wait, rendering (step 1/2)...');
parfor iDir = 1:nDirs
    hbar.iterate(1);
    BRIR_TimeDataTemp = Synthesize_SDM_Binaural(SRIR_data, BRIR_data, HRTF_TransL, HRTF_TransR, BRIR_data.Directions(iDir,:),0);
    BRIR_TimeDataTemp = ModifyReverbSlope(BRIR_data, BRIR_TimeDataTemp, OriginalT30, DesiredT30, RTFreqVector);
    BRIR_Early(:,:,iDir) = BRIR_TimeDataTemp;
end
close(hbar)

% Render late reverb
BRIR_full = Synthesize_SDM_Binaural(SRIR_data, BRIR_data, HRTF_TransL, HRTF_TransR, [0 0],1);
BRIR_full = ModifyReverbSlope(BRIR_data, BRIR_full, OriginalT30, DesiredT30, RTFreqVector);

% Remove leading zeros
[BRIR_Early, BRIR_full] = removeInitialDelay(BRIR_Early,BRIR_full,-20,BRIR_data.MixingTime*BRIR_data.fs);

% Split the BRIR
[early_BRIR, late_BRIR, DS_BRIR, ER_BRIR]  = split_BRIR(BRIR_Early, BRIR_full, BRIR_data.MixingTime, BRIR_data.fs, 256);

% -----------------------------------------------------------------------
% 5. Apply AP processing for the late reverb

% AllPass filtering for the late reverb (increasing diffuseness and
% smoothing out the EDC)
allpass_delays = [37 113 215 347];                      % in samples
allpass_RT = [0.1 0.1 0.1 0.1];                         % in seconds

for iAllPass=1:3
    late_BRIR(:,1) = allpass_filter(late_BRIR(:,1),allpass_delays(iAllPass) , [0.1], 48e3);
    late_BRIR(:,2) = allpass_filter(late_BRIR(:,2),allpass_delays(iAllPass) , [0.1], 48e3);
end

% -----------------------------------------------------------------------
% 6. Save BRIRs

hbar = parfor_progressbar(nDirs,'Please wait, saving (step 2/2)...');
parfor iDir = 1:nDirs
    SaveBRIR(SRIR_data, BRIR_data, DS_BRIR(:,:,iDir), early_BRIR(:,:,iDir), ER_BRIR(:,:,iDir), late_BRIR,BRIR_data.Directions(iDir,:));
    hbar.iterate(1);
end
SaveRenderingStructs(SRIR_data, BRIR_data)
close(hbar)



    