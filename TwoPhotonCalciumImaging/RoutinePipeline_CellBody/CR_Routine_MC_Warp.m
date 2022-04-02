%% Motion correction
% MC: Aki's code, Frontiers in Neuroinformatics, DOI:10.3389/fninf.2018.00098

clear all
close all
clc

tic
Rig = 'BScope1';
Animal = 'CR_4412583-LR';
Dates = {'210927','210928','210929','210930','211001','211002'};
Session = 'Rec_1'; % Recording session
signal_channels = 1; % red = 2,green = 1
isPiezo = false; % is alternating between two focal planes

for ii = 1:length(Dates)
    Date = Dates{ii};
    CR_Routine_MC_Daily(Rig,Date,Animal,Session,signal_channels,isPiezo)
end
toc

%% Warp
% PatchWarp if necessary, especially for axonal imaging
% https://doi.org/10.1101/2021.11.10.468164

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
% check which code I'm using
which ecc_old
which spatial_interp

Animals = {'CR_3575265-LR'};
Dates = {'190323','190324','190325','190326','190327','190328'};

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

