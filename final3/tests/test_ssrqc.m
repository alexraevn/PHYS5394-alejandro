%% Test the ssrqc function from glrtqcsig4pso.m

samplFreq = 1024;
interval = 1/samplFreq;
xVec = 0:interval:1;
nSamples = length(xVec);
% qc parameters
x = [2,0.3,10];

% psd handle
tPSD = @(f) (f>=10 & f<=200).*(f-10).*(200-f)/10000;

kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1));
psdVec = tPSD(posFreq);
% plot(freqVec,psdVec);

% generate the data

outNoise = statgaussnoisegen(nSampl,[posFreq(:),psdVec(:)],500,samplFreq);

qc = crcbgenqcsig(xVec,snr,x);

dataY = outNoise+qc;

% parameter structure
params = struct('dataX',xVec,...
                'dataXSq',xVec.^2,...
                'dataXCb',xVec.^3,...
                'dataY',dataY,...
                'samplFreq',samplFreq,...
                'psdVec',psdVec,...
                'snr',snr);

llr = glrtqcsig(dataY, samplFreq, psdVec, x);

ssr = ssrqc(x,params);

function ssrVal = ssrqc(x,params)
%Generate normalized quadratic chirp
phaseVec = x(1)*params.dataX + x(2)*params.dataXSq + x(3)*params.dataXCb;
qc = sin(2*pi*phaseVec);
% Inner product of signal with itself
normSigSqrd = innerprodpsd(qc,qc,params.samplFreq,params.psdVec);
% Normalization factor
normFac = params.snr/sqrt(normSigSqrd);
% Normalize qc signal
qc = normFac*qc;

%Compute fitness
ssrVal = -innerprodpsd(params.dataY,qc,params.samplFreq,params.psdVec)^2;
end