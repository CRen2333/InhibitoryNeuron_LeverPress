%% Motion correction and warp
clear all
close all
clc

% Motion correction
tic
Rig = 'BScope1';
Animal = 'CR_4412583-LR';
Dates = {'210927','210928','210929','210930','211001','211002'};
Session = 'Rec_1';
signal_channels = 1; % red = 2,green = 1
isPiezo = false;
for ii = 1:length(Dates)
    Date = Dates{ii};
    CR_Routine_MC_Daily(Rig,Date,Animal,Session,signal_channels,isPiezo)
end
toc

% Warp, patch
clear;
path_n = path;
if ~isempty(strfind(path_n,'Z:\People\Chi\warp for Chi'))
    rmpath(genpath('Z:\People\Chi\warp for Chi'));
end
path_n = path;
if isempty(strfind(path_n,'C:\Lab\Code\TwoP_IN\Warp'))
    disp('Add patch warping path...');
    addpath(genpath('C:\Lab\Code\TwoP_IN\Warp\Old'));
end
if ~isempty(strfind(path_n,'C:\Lab\Code\TwoP_IN\Warp\Old'))
    disp('Remove old warping path...')
    rmpath(genpath('C:\Lab\Code\TwoP_IN\Warp\Old'));
end
which ecc_old
which spatial_interp

Animals = {'CR_4412583-LR'};
Dates = {'210927','210928','210929','210930','211001','211002'};

for ii = 1:length(Animals)
    Animal = Animals{ii};
    for jj = 1:length(Dates)
        Date = Dates{jj};
        islocal = false;
        foldername = fullfile('F:\Data\MotionCorrection\',Animal,Date,'Rec_1\Z1\motioncorrected_tiff');
        local2savefoldername = fullfile('F:\Data\MotionCorrection\',[Animal '_Warp'],Date,'Rec_1\Z1\motioncorrected_tiff');
        rl_RemovingPix = 0;
        block_size = 6;
        overlap_npix = 10;
        n_split4warpinit = 1;
        
        CR_warp_Npatches_RH4tscc(islocal, foldername, local2savefoldername, rl_RemovingPix, block_size, overlap_npix, n_split4warpinit);
    end
end