
clear all
clc
Path='G:\WCJ-Data\双压杆社会抉择-WCJ\ACC数据\Figure 电生理\ACC电生理分析\DACC-plexondata\处理后的nex5文件\test-t\new\width_150ms_bin_50mss1\GPFA\';

file=dir(fullfile([Path '*mat']))
for ii=1:length(file)
    Name=file(ii).name
dat=load([Path Name]).datam;

% Results will be saved in mat_results/runXXX/, where XXX is runIdx.
% Use a new runIdx for each dataset.

runIdx =ii;

 method = 'gpfa';

% Select number of latent dimensions
xDim = 8;
% NOTE: The optimal dimensionality should be found using 
%       cross-validation (Section 2) below.

% If using a two-stage method ('fa', 'ppca', or 'pca'), select
% standard deviation (in msec) of Gaussian smoothing kernel.
kernSD = 30;
% NOTE: The optimal kernel width should be found using 
%       cross-validation (Section 2) below.
% Extract neural trajectories
result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD);
% NOTE: This function does most of the heavy lifting.
% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
% NOTE: The importance of orthnormalization is described on 
%       pp.621-622 of Yu et al., J Neurophysiol, 2009.
for i=1:length(dat)
    seqTrain(i).event=dat(i).event;
end

%plot3D(seqTrain, 'xsm', 'dimsToPlot', i);
% Plot neural trajectories in 3D space
plot3D(seqTrain, 'xorth', 'dimsToPlot', i);
savefig([Path Name(1:end-3) 'fig'])
save([Path Name(1:end-3) 'Seqtrian.mat'],'seqTrain')
end