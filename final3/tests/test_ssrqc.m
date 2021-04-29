%% Test the ssrqc function from glrtqcsig4pso.m

samplFreq = 1024;
nSamples = 2048;
% interval = 1/samplFreq;
xVec = (0:(nSamples-1))/sampFreq;
% nSamples = length(xVec);
% qc parameters
x = [2,0.3,10];
snr = 10;
% psd handle
tPSD = @(f) (f>=10 & f<=200).*(f-10).*(200-f).*((sin(f/6)*0.3+1))/10000 + 1;

kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdVec = tPSD(posFreq);
% plot(posFreq,psdVec);

% generate the data

outNoise = statgaussnoisegen(nSamples,[posFreq(:),psdVec(:)],500,samplFreq);

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

llr = glrtqcsig(dataY, samplFreq, psdVec, x); %Nan?
ssr = ssrqc(x,params);

% disp(llr);
% disp(ssr);

if round(llr,3) == -round(ssr,3)
    disp('glrtqcsig and ssrqc agree.')
end

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