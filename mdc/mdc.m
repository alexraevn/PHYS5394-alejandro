%% Significance estimation and parameter estimation of mock data
% parameter search range and length of ranges
a1 = [40,100]; a1l = a1(2)-a1(1);
a2 = [1,50];   a2l = a2(2)-a2(1);
a3 = [1,15];   a3l = a3(2)-a3(1);

% read TRAINIGDATA.MAT and ANALYSISDATA.MAT
trainingdatamat = load('TrainingData.mat');
dTrain = trainingdatamat.trainData;
sampFreq = trainingdatamat.sampFreq; % Hz
dataVec = load('analysisData.mat').dataVec;

% gen time vector
nSamples = length(dataVec);
dataLen = nSamples/sampFreq;
dataX = (0:(nSamples-1))/sampFreq;

% estimate PSD
%SDM************
%[pxx,posFreq] = pwelch(dTrain, 80,[],[],sampFreq);
[pxx,posFreq] = pwelch(dTrain, 1024,[],[],sampFreq);
%***************
% interpolate psd and freq vec to get length 1024
psdVec = interp1(1:length(pxx),pxx,linspace(1,length(pxx),1025),'cubic');
posFreq = interp1(1:length(posFreq),posFreq,linspace(1,length(posFreq),1025),'linear');

%% Set up matrix of parameters to get get GLRT values for
% Number of parameters to test (tn^3).
tn = 4;
% vectors of each parameter I'm testing of length tn
ta1 = a1(1):floor(a1l/(tn-1)):a1(2);
ta2 = a2(1):floor(a2l/(tn-1)):a2(2);
ta3 = a3(1):floor(a3l/(tn-1)):a3(2);
% vectors to use iteratively to populate MASTERPARS
stack3 = ones(tn,3);
stack2 = ones(tn^2,3);
% matrix with tn^3 qc sets of parameters to test
masterPars = ones(tn^3,3);

for q1 = 1:tn
    for q2 = 1:tn
        for q3 = 1:tn
              stack3(q3,:) = [ta1(q1),ta2(q2),ta3(q3)];
        end
        stack2((q2*tn-tn+1):q2*tn,:) = stack3;
    end
    masterPars((q1-1)*tn^2+1:q1*tn^2,:) = stack2;
end

%% Calculate the GLRT
% First I'm (t)esting 3 points in each qc parameter range
% % Lowest significance for 3^3 appeared around [70,49,8] at 0.2.
% % Lowest significance for 4^3 appeared around [80,49,1] at 0.18.
% % Then I searched only around [[65,85],[40,60],[1,15]];
% % % Lowest significance around [77,52,5] at 0.037.
%% I belive there is a qc sig in the data from such low significance.
%
% m noise realizations
m = 1000;
% glrtH0 vector of the GLRT values for m noise realizations, for tn^2 parameters
glrtH0 = zeros(tn^3,m);
allglrt = zeros(tn^3); % all glrt values
allsign = zeros(tn^3); % all significances
% Manually create fir2 filter once to save time
sqrtPSD = sqrt(psdVec);
b = fir2(100,posFreq/(sampFreq/2),sqrtPSD);

% uncomment to estimate glrt and significance of the paramenter ranges
for t = 1:(tn^3)
    for r = 1:m % m realizations
        % Generate WGN realization to pass through filter
        inNoise = randn(1,nSamples);
        noiseVec = sqrt(sampFreq)*fftfilt(b,inNoise);
        glrtH0(t,r) = glrtqcsig(noiseVec,sampFreq,psdVec,masterPars(t,:));
    end
    allglrt(t) = glrtqcsig(dataVec, sampFreq, psdVec, masterPars(t,:));
    allsign(t) = sum(glrtH0(t,:)>=allglrt(t))/m;
end

% I manually studied graph
% plot(allsign);
[mSign,mIndex] = min(allsign);
disp(['Lowest significance=',num2str(mSign(1)),newline,...
'Parameters a1=',num2str(masterPars(mIndex(1),1)),...
         '; a2=',num2str(masterPars(mIndex(1),2)),...
         '; a3=',num2str(masterPars(mIndex(1),3))]);

%% Call PSO
% parameter structrue for glrtqc4pso
nRuns = 8;
nSteps = 2000;
snr = 1;
inParams = struct('dataX',dataX,...
                  'dataXSq',dataX.^2,...
                  'dataXCb',dataX.^3,...
                  'dataY',dataVec,...
                  'samplFreq',sampFreq,...
                  'psdVec',psdVec,...
                  'snr',snr,...
                  'rmin',[a1(1),a2(1),a3(1)],...
                  'rmax',[a1(2),a2(2),a3(2)]);

psoParams = struct('steps',nSteps);
% uncomment to run PSO
outStruct = glrtqcpso(inParams,psoParams,nRuns);

%% Estimation of SNR
% Norm of best signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(outStruct.bestSig,outStruct.bestSig,sampFreq,psdVec);
% Normalize best signal to specified SNR
sigVec = snr*outStruct.bestSig/sqrt(normSigSqrd);

% Obtain LLR values under H0 for bestSig
flatness = 100; % time samples to truncate at the beginning of noise
nH0Data = 500;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    inNoise = randn(1,nSamples+flatness);
    noiseVec = sqrt(sampFreq)*fftfilt(b,inNoise);
    noiseVec = noiseVec(flatness+1:end);
    llrH0(lp) = innerprodpsd(noiseVec,sigVec,sampFreq,psdVec);
end
%Obtain LLR for multiple data (=signal+noise) realizations
% nH1Data = 500;
% llrH1 = zeros(1,nH1Data);
% for lp = 1:nH0Data
%     inNoise = randn(1,nSamples+flatness);
%     noiseVec = sqrt(sampFreq)*fftfilt(b,inNoise);
%     noiseVec = noiseVec(flatness+1:end);
%     % Add normalized signal
%     H1data = noiseVec + sigVec;
%     llrH1(lp) = innerprodpsd(H1data,sigVec,sampFreq,psdVec);
% end
llrH1 = innerprodpsd(dataVec,sigVec,sampFreq,psdVec);

%Signal to noise ratio estimate
estSNR = (llrH1-mean(llrH0))/std(llrH0);
bestSig = estSNR*sigVec;

%% Presentation

figure;
hold on;
plot(dataX,dataVec,'Color','Magenta');
% plot(dataX,sigVec,'Color','Red','LineWidth',2.0);
plot(dataX,bestSig,'LineStyle','-','LineWidth',1.0, 'Color',[76,153,0]/255);
title('Data and Estimated Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Data','Estimated Signal');
disp(['PSO estimated parameters: a1=',num2str(outStruct.bestQcCoefs(1)),...
                              '; a2=',num2str(outStruct.bestQcCoefs(2)),...
                              '; a3=',num2str(outStruct.bestQcCoefs(3)),...
                              newline,'Estimated SNR=',num2str(estSNR)]);

% %% Spectrograms
% % Declare length of windows and their overlap
% winLen = 0.08;
% ovrlp  = 0.07;
% 
% % Convert to integer number of samples
% winLenSamp = floor(winLen*sampFreq);
% ovrlpSamp  = floor(ovrlp*sampFreq);
% 
% % Generate spectrogram of testData.txt
% [Sin,Fin,Tin] = spectrogram(dataVec, winLenSamp, ovrlpSamp, [], sampFreq);
% 
% % Plot data spectrogram
% figure;
% imagesc(Tin, Fin, abs(Sin)); axis xy;
% title('Spectrogram of Data');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% colorbar;







