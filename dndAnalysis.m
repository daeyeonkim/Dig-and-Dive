%% Analysis of DnD behavior %%
% Develped by Daeyeon Kim    %
% Dec. 11st, 2013            %


clear
close all

%% SET FOLDERS --------------------------------------------------------- %%
% input directory
DirHome = '/Users/daeyeonkim/Documents/Box Sync/MATLAB/projects/DnD/codes/';
DirImport = '/Users/daeyeonkim/Documents/Box Sync/MATLAB/projects/DnD/data/sample/';

DirExport = strcat(DirImport,'analyzed/');
mkdir(DirExport);

cd(DirImport);

% array folders
FolderArray = dir(strcat(DirImport,'tracking*'));
% chose experimental conditions
conditions = [1];
maxFolder = length(conditions);

%% SET PARAMETERS ------------------------------------------------------ %%
% set constant variables ------------------------------------------------ %
% bin depth of larva and time
zbins = 0.1;  % best: 0.1
tbins = 4; % best: 3 or 4
% thresholds for digging and diving motions
thDig = -0.0;
thDive = [-3.6 -3.6];

maxLarvae = 12; % maximum number of larvae
analTime = 4; % assay time in min
binSliding = 60/60; % in min
fps = 2; % frame per sec (Hz)

% figure options
fontSize = 10;
rotation = 30;
RasterON = 1; % ON=1, OFF=0
EthogramON = 0; % ON=1, OFF=0
BoxplotON = 1;
SaveDataON = 1;
SaveStatON = 0;

% define conditions
geno = 'CantonS-olfaction';

cond{1}='agar01';
%cond{2}='EtA+';

% define matrixes
quan{1} = 'TotalDivingTime';
quan{2} = 'TotalDiggingTime';
quan{3} = 'TotalSurfacingTime';
quan{4} = 'Risk-taking index';

% define variables ------------------------------------------------------ %
AlltotTdive = [];
AlltotTdig = [];
AlltotTsurf = [];
AllSingleDiveT = [];
AllMaxDiveDepth = [];
AllRI = []; % All risk-taking index

AllnDive = []; % number of dives
AllfDive = []; % freqeuncy of dives
meanFdive = []; % mean dive frequency
semFdive = []; % sem of dive frequecy
meanIntDive = []; % mean time interval between dives
semIntDive = []; % sem of time interval

meanVelDiveAll = [];
semVelDiveAll = [];
workDiveAll = [];
diveSpeedAll = [];
meanVelSurfAll = [];
semVelSurfAll = [];

groupLarvae = [];
groupDive = [];

pDiveDigAll = [];
pDigDiveAll = [];
pDigSurfAll = [];
pSurfDigAll = [];
pDiveDiveAll = [];
pDigDigAll = [];
pSurfSurfAll = [];
pDiving = [];
pDigging = [];
pSurfacing = [];

trMatrix = {};

Labels = {};
nLarvae = [];

cc = {'black','red','blue','green','magenta','cyan','yellow','white','green','red','blue','blue'};

%% loading and sorting data -------------------------------------- %%
for folder = 1:maxFolder
    % set the current folder
    cd(DirImport);
    FolderName = FolderArray(conditions(folder)).name;
    cd(FolderName);
    load('trackingDataAll.mat'); % load data file
    
    % define constant variables ----------------------------------------- %
    IdxLarva = 1:length(trackingDataAll);
    nLarvae(folder) = IdxLarva(end);
    Labels{1,folder} = cond{folder};
    
    % define variables -------------------------------------------------- %
    maxs=[];
    mins=[];
    idxmaxs=[];
    idxmins=[];
        
    totTdive = [];
    totTdig = [];
    totTsurf = [];
    
    nDive = [];
    intDive = [];
   
    sumDive = [];
    sumDig = [];
    sumRest = [];
    allDive = [];
    allIniDive = [];
    allDigAfDive = [];
    allSurfAfDive = [];
    allMdive = [];
    allIntDive = [];
    allBeforeDiveInterval = [];
    
    DiveProb = [];
    DigProb = [];
    SurfProb = [];
    idxDiveMinAll = [];
    idxDigMinAll = [];
    idxSurfMinAll = [];
    
    allDig = [];
    allSurf = [];
    
    maxDiveDepth = [];
    meanVelDive = [];
    workDive = [];
    meanVelSurf = [];

    anotModeAll = [];
    
    perDive = [];
    perDig = [];
    perRest = [];
    AllPerDive = [];
    AllPerDig = [];
    AllPerRest = [];
    
    IntegTotal = [];
    IntegDiveSum = [];
    RI = [];  % risk-taking index
    PI = []; % preference index
    
    % correction for depth of larva by replacing NaN with previous depth,
    % which is not NaN ------------------------------------------------- %
    for i=IdxLarva
        depth = trackingDataAll{i}.verticalCOM(1:analTime*60*fps); % define the depth of larva
        depth0 = length(depth);
        if (depth0<analTime*60*fps)
            depth(depth0+1:analTime*60*fps) = depth(end);
        end
        for j = 1:length(depth)
            if (isnan(depth(j))==1)
                if (j==1) % if NaN is first, take the first depth of larva
                   depth(j) = depth(min(find(isnan(depth)==0)));
                else
                   depth(j)=depth(j-1); % replace NaN with the previous #
                end
            end
        end
        totTdive(i) = length(find(depth <= thDive(folder)))/fps; % total dive time
    end
    
    % set variables ---------------------------------------------------- %
    time = (1:1:length(depth))/60/fps; % assay time, unit of min
    totTime = time(end);
    % number of transitions
    nDiveDig = zeros(6,maxLarvae); % transitions Diving(1) --> Digging(2)
    nDigDive = zeros(6,maxLarvae); % transitions Digging(2) --> Diving(1)
    nDigSurf = zeros(6,maxLarvae); % transitions Digging(2) --> Surfacing(3)
    nSurfDig = zeros(6,maxLarvae); % transitions Surfacing(3) --> Digging(2)
    nDiveDive = zeros(1,maxLarvae); % 1 --> 1
    nDigDig = zeros(1,maxLarvae); % 2 --> 2
    nSurfSurf = zeros(1,maxLarvae); % 3 --> 3
    % prob of transitions
    pDiveDig = zeros(6,maxLarvae);
    pDigDive = zeros(6,maxLarvae);
    pDigSurf = zeros(6,maxLarvae);
    pSurfDig = zeros(6,maxLarvae);
    pDiveDive = zeros(1,maxLarvae);
    pDigDig = zeros(1,maxLarvae);
    pSurfSurf = zeros(1,maxLarvae);
    
    % sort larvae by dive time ----------------------------------------- %
    
    [s1, s2] = sort(totTdive,'descend'); % sort dive time in descend
    IdxLarva = s2; % re-assign index of larva using sorted larvae
    n=0; % local index of larva
    % correction of depth again ----------------------------------------- %
    for i=IdxLarva
        n = n+1;
        depth = trackingDataAll{i}.verticalCOM(1:analTime*60*fps);
        x = trackingDataAll{i}.horizontalCOM(1:analTime*60*fps);
        depth0 = length(depth);
        if (depth0<analTime*60*fps)
            depth(depth0+1:analTime*60*fps) = depth(end);
            x(depth0+1:analTime*60*fps) = x(end);
        end
        
        for j = 1:length(depth)
            if (isnan(depth(j))==1)
                if (j==1)
                   depth(j)= depth(min(find(isnan(depth)==0)));
                else
                depth(j) = depth(j-1);
                end
            end
            if (isnan(x(j))==1)
                if (j==1)
                    x(j)= x(min(find(isnan(x)==0)));
                else
                    x(j)=x(j-1);
                end
            end
        end
           
    % calculate speed -------------------------------------------------- %
    dz = depth(2:end)-depth(1:end-1);
    dx = x(2:end)-x(1:end-1);
    vel(i,:) = [sqrt(dx.^2 + dz.^2)/fps , sqrt(dx(end)^2 + dz(end)^2)/fps];
    % add the same speed value to the end to get the same length of vector with depth       
        
    %% ANOTATE EACH MODE ----------------------------------------------- %%         
    anotMode = []; % anotation of each mode
    idxDive = find(depth <= thDive(folder)); % index in diving mode
    idxDig = find(depth <= thDig & depth > thDive(folder));  % index in digging mode
    idxSurf = find(depth  > thDig); % index in surfacing mode
    
    anotMode(idxDive) = 1; % anotation of diving mode as '1'
    anotMode(idxDig) = 2; % anotation of digging mode as '2'
    anotMode(idxSurf) = 3; % anotation of surfacing mode as '3'
    anotModeAll = cat(1,anotModeAll,anotMode);
    
    %% CALCULATE ETHOGRAM ---------------------------------------------- %%
    % set variables ----------------------------------------------------- %
    trans = []; % whole transitions
    nTrans = 0;
    tSelf = []; % transitions itself
    tModes = []; % transitions between modes
   
    % count transitions between modes ----------------------------------- %    
    trans = diff(anotMode); % define transitions
    nTrans = length(trans); % number of transitions
    tSelf = find(trans == 0);
    % self transitions
    nDiveDive(n) = length(find(anotMode(tSelf) == 1));
    nDigDig(n) = length(find(anotMode(tSelf) == 2));
    nSurfSurf(n) = length(find(anotMode(tSelf) == 3));
    tModes = find(trans ~= 0);
    % transitions between different modes
    for subT=0:2
    tModeSub = tModes(tModes > subT*5*60 & tModes <= (subT+1)*5*60);
    for k = 1:length(tModeSub)
        if (trans(tModeSub(k)) > 0 && anotMode(tModeSub(k)) == 1)
            nDiveDig(subT+1,n) = nDiveDig(subT+1,n) + 1; % increase # transition 1-2
        elseif (trans(tModeSub(k)) > 0 && anotMode(tModeSub(k)) == 2)
            nDigSurf(subT+1,n) = nDigSurf(subT+1,n) + 1; % increase # transition 2-3
        elseif (trans(tModeSub(k)) < 0 && anotMode(tModeSub(k)) == 2)
            nDigDive(subT+1,n) = nDigDive(subT+1,n) + 1; % increase # transition 2-1
        elseif (trans(tModeSub(k)) < 0 && anotMode(tModeSub(k)) == 3)
            nSurfDig(subT+1,n) = nSurfDig(subT+1,n) + 1; % increase # transition 3-2
        else
        end
    end
    % prob of transitions (only counting tModes)
    pDiveDive(subT+1,n) = 1;
    pDigDig(subT+1,n) = 1;
    pSurfSurf(subT+1,n) = 1;
    
 
    pDiveDig(subT+1,n) = nDiveDig(subT+1,n)/length(tModeSub);
    pDigDive(subT+1,n) = nDigDive(subT+1,n)/length(tModeSub);
    pDigSurf(subT+1,n) = nDigSurf(subT+1,n)/length(tModeSub);
    pSurfDig(subT+1,n) = nSurfDig(subT+1,n)/length(tModeSub);
    
    pDiveDig(isnan(pDiveDig)==1) = 0;
    pDigDive(isnan(pDigDive)==1) = 0;
    pDigSurf(isnan(pDigSurf)==1) = 0;
    pSurfDig(isnan(pSurfDig)==1) = 0;
    
    end
    
% prob of transitions (conuting all transitions)
%     pDiveDive(n) = nDiveDive(n)/nTrans;
%     pDigDig(n) = nDigDig(n)/nTrans;
%     pSurfSurf(n) = nSurfSurf(n)/nTrans;
%     
%     pDiveDig(n) = nDiveDig(n)/nTrans;
%     pDigDive(n) = nDigDive(n)/nTrans;
%     pDigSurf(n) = nDigSurf(n)/nTrans;
%     pSurfDig(n) = nSurfDig(n)/nTrans;
    

    %% CALCULATE TIME SPENT IN EACH MODE ------------------------------- %%
    % set variables ----------------------------------------------------- %
    trDive = []; % transitions in derivative of idxDive
    trDig = []; % transitions in derivative of idxDig
    trSurf = []; % transitions derivative of idxSurf
    % variables for start and end of each mode
    idxStartDive = [];
    idxEndDive = [];
    idxStartDig = [];
    idxEndDig = [];
    idxStartSurf = [];
    idxEndSurf = [];
    % find start and end of each mode ----------------------------------- %
    % find start and end of diving time
    if (isempty(idxDive)==1) % in case of no dives
        %idxDive = 0;
        idxDive = nan;
    else
        trDive = find(diff(idxDive)>1)+1;
        if (isempty(trDive)==1) % in case all modes are divings
           idxStartDive = idxDive(1);
           idxEndDive = idxDive(end);
        else % in case of descrete dives
           idxStartDive = cat(2,idxDive(1),idxDive(trDive));
           idxEndDive = cat(2,idxDive(trDive-1),idxDive(end));
        end
    end
    % find start and end of diggning time
    if (isempty(idxDig)==1)
        %idxDig = 0;
        idxDig = nan;
    else
        trDig = find(diff(idxDig)>1)+1;
        if (isempty(trDig)==1)
           idxStartDig = idxDig(1);
           idxEndDig = idxDig(end);
        else
           idxStartDig = cat(2,[idxDig(1)],idxDig(trDig));
           idxEndDig = cat(2,idxDig(trDig-1),idxDig(end));
        end
    end    
    % find start and end of surfacing time
    if (isempty(idxSurf)==1)
        %idxSurf = 0;
        idxSurf = nan;
    else
        trSurf = find(diff(idxSurf)>1)+1;
        if (isempty(trSurf)==1)
           idxStartSurf = idxSurf(1);
           idxEndSurf = idxSurf(end);
        else
           idxStartSurf = cat(2,idxSurf(1),idxSurf(trSurf));
           idxEndSurf = cat(2,idxSurf(trSurf-1),idxSurf(end));
        end
    end
    
    % calculate single time for each mode ------------------------------- %   
    tDive = idxEndDive - idxStartDive + 1; % including the first frame
    tDig = idxEndDig - idxStartDig + 1;
    tSurf = idxEndSurf - idxStartSurf + 1;
       
    % discard modes having no time durations ---------------------------- %
    idxNull = tDive==0;
    tDive(idxNull)=[];
    idxStartDive(idxNull)=[];
    idxEndDive(idxNull)=[];
    
    idxNull = tDig==0;
    tDig(idxNull)=[];
    idxStartDig(idxNull)=[];
    idxEndDig(idxNull)=[];
    
    idxNull = tSurf==0;
    tSurf(idxNull)=[];
    idxStartSurf(idxNull)=[];
    idxEndSurf(idxNull)=[];
    
    % calculate total time for each mode -------------------------------- %
    totTdive(i) = sum(tDive);
    totTdig(i) = sum(tDig);
    totTsurf(i) = sum(tSurf);
    
    % calculate time interval between dives ----------------------------- %
    intDive = idxStartDive(2:end) - idxEndDive(1:end-1);
    %intDive = idxStartDive(2:end) - idxStartDive(1:end-1);
    
        %% CALCULATE RISK-TAKING INDEX ------------------------------------- %
    % Calculate total integral z(t) over t

    IntegDive = [];
    IntegSurf = [];
    IntegDig = [];
    IntegTotal(i) = trapz(depth)-thDive(folder)*time(end)*60*fps;
    %IntegTotal(i) = trapz(depth);
    % Calculate dive integral z(t)@dive over t
    
    for j=1:length(idxStartDive)
        if (idxStartDive(j) == idxEndDive(j))
            IntegDive(j) = 0;
        else
            IntegRange = idxStartDive(j):1:idxEndDive(j);
            IntegDive(j) = trapz(IntegRange,depth(IntegRange));
        end
    end
    
    for j=1:length(idxStartDig)
        if (idxStartDig(j) == idxEndDig(j))
            IntegDig(j) = 0;
        else
            IntegRange = idxStartDig(j):1:idxEndDig(j);
            IntegDig(j) = trapz(IntegRange,depth(IntegRange));
        end
    end
   
%     for j=1:length(idxStartSurf)
%         if (idxStartSurf(j) == idxEndSurf(j))
%             IntegSurf(j) = 0;
%         else
%             IntegRange = idxStartSurf(j):1:idxEndSurf(j);
%             IntegSurf(j) = trapz(IntegRange,depth(IntegRange));
%         end
%     end
    
    IntegDiveSum(i) = sum(IntegDive)-thDive(folder)*totTdive(i);
    %IntegSurfSum(i) = sum(IntegSurf);
    IntegDigSum(i) = (0-thDive(folder))*totTdig(i);
    IntegSurfSum(i) = 2*totTsurf(i);
    
    % Risk-taking index
    %RI(i) = IntegDiveSum(i)/(IntegTotal(i)-IntegDiveSum(i));
    %RI(i) = IntegDiveSum(i)/(IntegTotal(i));
    
    %RI(i) = IntegDiveSum(i)/(IntegTotal(i)-2*IntegDiveSum(i));
    
    RI(i) = -IntegDiveSum(i)/-sum(IntegDig);  % recent one
    
    %RI(i) = -IntegDiveSum(i)/(time(end)*60*fps*(12+thDive(folder)));  % recent one
    %RI(i) = (IntegTotal(i)+IntegDiveSum(i))/(IntegTotal(i)-1*IntegDiveSum(i));  % recent one
    %RI(i) = IntegDiveSum(i)/(IntegTotal(i)-1*IntegDiveSum(i));
    %RI(i) = IntegTotal(i);
    %RI(i) = (-IntegDiveSum(i)/(time(end)*60*fps*(12+thDive(folder))) + IntegSurfSum(i)/(time(end)*60*fps*2))/2;  % recent one
    %RI(i) = (-IntegDiveSum(i) + IntegSurfSum(i))/IntegDigSum(i);  % recent one
    %RI(i) = (-IntegDiveSum(i) + IntegSurfSum(i))/(time(end)*60*fps*(14+thDive(folder)));  % recent one
    %% ----
    PI(i) = (totTsurf(i) - totTdive(i))/(totTsurf(i)+totTdive(i)); % preference index  
    
    % Get max. dive depth ----------------------------------------------- %
    for j=1:length(tDive)
    maxDiveDepth = cat(2,maxDiveDepth,-min(depth(idxStartDive(j):idxEndDive(j))));
    end
        
    % calculate number of dives ---------------------------------------- %
    nDive(n) = length(tDive);
    if (length(tDive)>1)
       allMdive = cat(2,allMdive,tDive(1:end-1));
    end

    % calculate dive speed --------------------------------------------- %
    if (isempty(tDive)==1)
    else
       for xx=1:length(tDive)
           idxlocal = idxStartDive(xx):idxEndDive(xx);
            velDive = nanmean(vel(i,idxlocal));
            %velDive(velDive>1)=[]; % get rid of very fast speed (wrong!)
          if (tDive(xx)==1)
             wDive = velDive^2;
          else
             wDive = trapz(idxlocal,vel(i,idxlocal).^2);
          end  
             meanVelDive = cat(2,meanVelDive,velDive);
             workDive = cat(2,workDive,wDive);
        end
    end
    
    % calculate surf speed ---------------------------------------------- %
    if (isempty(tSurf)==1)
    else
       for xx=1:length(tSurf)
           idxlocal = idxStartSurf(xx):idxEndSurf(xx);
            velSurf = nanmean(vel(i,idxlocal));
            %velDive(velDive>1)=[]; % get rid of very fast speed (wrong!)
%           if (tSurf(xx)==1)
%              wDive = velDive^2;
%           else
%              wDive = trapz(idxlocal,vel(i,idxlocal).^2);
%           end  
             meanVelSurf = cat(2,meanVelSurf,velSurf);
             %workDive = cat(2,workDive,wDive);
        end
    end
      
%     idxStartDive(meanVelDive>1)=[];
%     idxEndDive(meanVelDive>1)=[];
%     tDive(meanVelDive>1)=[];
%     workDive(meanVelDive>1)=[];
%     meanVelDive(meanVelDive>1)=[];
  
    %% CALCULATE PROB. DIVE, DIGGING, AND SURFACING OVER TIME ---------------------------------- %%
%     anotOnlyDive = zeros(1,length(anotMode));
%     anotOnlyDive(find(anotMode == 1)) = 1;
%     probDive = anotOnlyDive/sum(anotOnlyDive);
    tbinmin = 1/fps/60;   % for all dive prob.:0.2, for meanDiveProb. : 0.3
%     binSliding = 1;
    idxDiveMin = idxDive/60/fps;
    idxDiveMinAll = cat(2,idxDiveMinAll,idxDiveMin);
    DiveRange = tbinmin:tbinmin:time(end);
    %DiveRange = 0:tbinmin:time(end);
%     DiveHist = [];
%     for j=1:length(time)
%         DiveHist(j) = numel(find(idxDiveMin >= time(j)/60 & idxDiveMin < time(j)/60+binSliding));
%     end
%     numDiveIdx = ones(1,length(DiveHist))*sum(DiveHist);
%     DiveProb(i,:) = DiveHist./numDiveIdx;    
    DiveProb(i,:) = histc(idxDiveMin,DiveRange)/length(idxDiveMin); 
    
    % for digging
    idxDigMin = idxDig/60/fps;
    idxDigMinAll = cat(2,idxDigMinAll,idxDigMin);
    DiveRange = 0:tbinmin:time(end);
    DigProb(i,:) = histc(idxDigMin,DiveRange)/length(idxDigMin); 
    
    % for surfacing
    idxSurfMin = idxSurf/60/fps;
    idxSurfMinAll = cat(2,idxSurfMinAll,idxSurfMin);
    DiveRange = 0:tbinmin:time(end);
    SurfProb(i,:) = histc(idxSurfMin,DiveRange)/length(idxSurfMin); 
    
    %% QUANTIFICATION OF DIGGING AND SURFACING TIME AFTER DIVE --------- %%

    % concatenating all sigle dives, digs, surfaces --------------------- %
    idxFirst = [];
    allDive = cat(2,allDive,tDive);
    allIniDive = cat(2,allIniDive,idxStartDive);
    allIntDive = cat(2,allIntDive,intDive);
    allDig = cat(2,allDig,tDig);
    allSurf = cat(2,allSurf,tSurf);
    
    diveTbefore = tDive(1:length(intDive));
    allBeforeDiveInterval = cat(2,allBeforeDiveInterval,[diveTbefore; intDive]);
    
   if (isempty(tDive)==0)
    
      for k=1:length(idxEndDive)
        if (k==length(idxEndDive))
%             idxFirst = find(idxStartDig >= idxEndDive(k),1);
%             if (isempty(idxFirst)==0)
%                 allDigAfDive = cat(2,allDigAfDive,tDig(idxFirst));
%             else
%                 allDigAfDive = cat(2,allDigAfDive,0);
%             end
%             %idxFirst = find(idxStartSurf >= idxEndDive(k));
%             idxFirst = find(idxStartSurf >= idxEndDive(k),1);
%             if (isempty(idxFirst)==0)
%                 %allSurfAfDive = cat(2,allSurfAfDive,sum(tSurf(idxFirst)));
%                 allSurfAfDive = cat(2,allSurfAfDive,tSurf(idxFirst));
%             else
%                 allSurfAfDive = cat(2,allSurfAfDive,0);
%             end
        else
            idxFirst = find(idxStartDig >= idxEndDive(k) & idxStartDig < idxStartDive(k+1),1);
            if (isempty(idxFirst)==0)
                allDigAfDive = cat(2,allDigAfDive,tDig(idxFirst));
            else
                allDigAfDive = cat(2,allDigAfDive,0);
            end
        
            %idxFirst = find(idxStartSurf >= idxEndDive(k) & idxStartSurf < idxStartDive(k+1));
            idxFirst = find(idxStartSurf >= idxEndDive(k) & idxStartSurf < idxStartDive(k+1),1);
            if (isempty(idxFirst)==0)
                %allSurfAfDive = cat(2,allSurfAfDive,sum(tSurf(idxFirst)));
                allSurfAfDive = cat(2,allSurfAfDive,tSurf(idxFirst));
            else
                allSurfAfDive = cat(2,allSurfAfDive,0);
            end
        end
      end
   end
   
   
    if (RasterON == 1)
    %% RASTER PLOTS --------------------------------------------------- %%
    % set variables ----------------------------------------------------- %
    yBarBottom = 0:nLarvae(folder)-1; % y-position of each rator plot (bottom)
    yBarTop = 1:nLarvae(folder); % y-position of each rator plot (top)
    yBar(n,:) = [yBarBottom(n) yBarTop(n) yBarTop(n) yBarBottom(n)];
    % ------------------------------------------------------------------- %
    
    fig1=figure(1);
   
    sp=subplot(3,2,folder);
    title(cond{folder})
    hold on 
    xlabel('Time (min)');
    ylabel('Larvae');
    %xlim([0 30])
%     set(sp,'XTick',[0:5:30])
%     set(sp,'XTickLabel',['0';'5';'10';'15';'20';'30'])
    ylim([0 nLarvae(folder)])
    
    % plot dives
    for j = 1:length(idxStartDive)
        xDive = [idxStartDive(j)-1 idxStartDive(j)-1 idxEndDive(j) idxEndDive(j)]/60/fps;
        h1 = fill(xDive,yBar(n,:),'g','EdgeColor','none');
    end
    % plot digs
    for j=1:length(idxStartDig)
        xDig = [idxStartDig(j)-1 idxStartDig(j)-1 idxEndDig(j) idxEndDig(j)]/60/fps;
        h2=fill(xDig,yBar(n,:),'b','EdgeColor','none');
    end
    % plot surfacings
    for j=1:length(idxStartSurf)
        xSurf = [idxStartSurf(j)-1 idxStartSurf(j)-1 idxEndSurf(j) idxEndSurf(j)]/60/fps;
        h3=fill(xSurf,yBar(n,:),'r','EdgeColor','none');
    end
    % plot lines between raters
    line([0,time(end)],[0,0],'color','k')
    line([0,time(end)],[n,n],'color','k')
    end
  end
    
    
    % save the results of ethogram -------------------------------------- %
    % combine prob of transitions between modes in all conditions
    pDiveDigAll = cat(1,pDiveDigAll,pDiveDig);
    pDigDiveAll = cat(1,pDigDiveAll,pDigDive);
    pDigSurfAll = cat(1,pDigSurfAll,pDigSurf);
    pSurfDigAll = cat(1,pSurfDigAll,pSurfDig);
    pDiveDiveAll = cat(1,pDiveDiveAll,pDiveDive);
    pDigDigAll = cat(1,pDigDigAll,pDigDig);
    pSurfSurfAll = cat(1,pSurfSurfAll,pSurfSurf);

    
    %% STATISTICS ---------------------------------------------------------- %%
    % mean transition matrix -----------------------------------------------  %
    trMatrixAll = [];
    for j=1:3
     
    trMatrixTemp(1,1) = nanmean(pDiveDive(j,1:nLarvae(folder)));
    trMatrixTemp(1,2) = nanmean(pDiveDig(j,1:nLarvae(folder)));
    trMatrixTemp(1,3) = 0; % no transitions from 1 --> 3
    trMatrixTemp(2,1) = nanmean(pDigDive(j,1:nLarvae(folder)));
    trMatrixTemp(2,2) = nanmean(pDigDig(j,1:nLarvae(folder)));
    trMatrixTemp(2,3) = nanmean(pDigSurf(j,1:nLarvae(folder)));
    trMatrixTemp(3,1) = 0; % no transitions from 3 --> 1
    trMatrixTemp(3,2) = nanmean(pSurfDig(j,1:nLarvae(folder)));
    trMatrixTemp(3,3) = nanmean(pSurfSurf(j,1:nLarvae(folder)));
    
    trMatrixAll = cat(2,trMatrixAll,trMatrixTemp);
    
    end
    
    trMatrix{folder} = trMatrixAll;
    
%     trMatrix{folder}(1,1) = mean(pDiveDiveAll(folder,1:nLarvae(folder)));
%     trMatrix{folder}(1,2) = mean(pDiveDigAll(folder,1:nLarvae(folder)));
%     trMatrix{folder}(1,3) = 0; % no transitions from 1 --> 3
%     trMatrix{folder}(2,1) = mean(pDigDiveAll(folder,1:nLarvae(folder)));
%     trMatrix{folder}(2,2) = mean(pDigDigAll(folder,1:nLarvae(folder)));
%     trMatrix{folder}(2,3) = mean(pDigSurfAll(folder,1:nLarvae(folder)));
%     trMatrix{folder}(3,1) = 0; % no transitions from 3 --> 1
%     trMatrix{folder}(3,2) = mean(pSurfDigAll(folder,1:nLarvae(folder)));
%     trMatrix{folder}(3,3) = mean(pSurfSurfAll(folder,1:nLarvae(folder)));
     
    % Group conditions ------------------------------------------------- %
    groupLarvae = cat(1,groupLarvae,ones(nLarvae(folder),1)*folder);
    groupDive = cat(1,groupDive,ones(length(allDive),1)*folder);
    
    %% CALCULATE DIVE FREQUENCY ---------------------------------------- %%
    AllnDive = cat(1,AllnDive,nDive');
    AllfDive = AllnDive/max(time); % all dive frequencies
    meanFdive(folder) = nanmean(AllfDive(groupLarvae==folder)); % mean dive frequency
    semFdive(folder) = nanstd(AllfDive(groupLarvae==folder))/sqrt(nLarvae(folder)); % sem
    
    %% CALCULATE TIME INTERVAL BETWEEN DIVES --------------------------- %%
    if (folder==1)
        nodor = allIntDive;
    else
        odor = allIntDive;
    end
    
    meanIntDive(folder) = nanmean(allIntDive);
    semIntDive(folder) = nanstd(allIntDive)/sqrt(length(allIntDive));
    
    %% GETTING All TOTAL TIME IN EACH MODE ---------------------------- %%
    AlltotTdive = cat(1,AlltotTdive,(totTdive/60/totTime*100)'); % all % time in dive mode
    AlltotTdig = cat(1,AlltotTdig,(totTdig/60/totTime*100)');   % all % time in dig mode
    AlltotTsurf = cat(1,AlltotTsurf,(totTsurf/60/totTime*100)');    % all % time in surf mode

    %% CALULATE DIVE, DIG, AND SURFACING PROB CROSS LARVAE OVER TIME ---------------------- %%
    MeanDiveProb(folder,:) = mean(DiveProb,1);
    SemDiveProb(folder,:) = std(DiveProb,1)/sqrt(size(DiveProb,1));
    
    DiveRange = tbinmin:tbinmin:time(end);
    DiveHist = [];
    for j=1:length(DiveRange)
        DiveHist(j) = numel(find(idxDiveMinAll >= DiveRange(j) & idxDiveMinAll < DiveRange(j)+binSliding));
    end
    numDiveIdx = ones(1,length(DiveHist))*nansum(DiveHist);
    DiveProbAll(folder,:) = DiveHist./numDiveIdx;     
    
%     DiveProbAll(folder,:) = histc(idxDiveMinAll,DiveRange)/length(idxDiveMinAll);
    
    MeanDigProb(folder,:) = mean(DigProb,1);
    SemDigProb(folder,:) = std(DigProb,1)/sqrt(size(DigProb,1));
    DigProbAll(folder,:) = histc(idxDigMinAll,DiveRange)/length(idxDigMinAll);
    
    MeanSurfProb(folder,:) = mean(SurfProb,1);
    SemSurfProb(folder,:) = std(SurfProb,1)/sqrt(size(SurfProb,1));
    SurfProbAll(folder,:) = histc(idxSurfMinAll,DiveRange)/length(idxSurfMinAll);
    
    %% Calculate prob. diving cross larvae
    for j=1:length(time)
        pDiving(folder,j)= length(find(anotModeAll(:,j)==1))/nLarvae(folder);
        pDigging(folder,j) = length(find(anotModeAll(:,j)==2))/nLarvae(folder);
        pSurfacing(folder,j) = length(find(anotModeAll(:,j)==3))/nLarvae(folder);
    end
    
    % applying moving average
    pDiving(folder,:) = tsmovavg(pDiving(folder,:),'s',binSliding*60*fps);
    pDiving(folder,binSliding/2*60*fps:time(end)*60*fps-binSliding/2*60*fps) = pDiving(folder,binSliding*60*fps:time(end)*60*fps);
    pDiving(folder,time(end)*60*fps-binSliding/2*60*fps:time(end)*60*fps) = nan;
    
    pDigging(folder,:) = tsmovavg(pDigging(folder,:),'s',binSliding*60*fps);
    pDigging(folder,binSliding/2*60*fps:time(end)*60*fps-binSliding/2*60*fps) = pDigging(folder,binSliding*60*fps:time(end)*60*fps);
    pDigging(folder,time(end)*60*fps-binSliding/2*60*fps:time(end)*60*fps) = nan;
    
    pSurfacing(folder,:) = tsmovavg(pSurfacing(folder,:),'s',binSliding*60*fps);
    pSurfacing(folder,binSliding/2*60*fps:time(end)*60*fps-binSliding/2*60*fps) = pSurfacing(folder,binSliding*60*fps:time(end)*60*fps);
    pSurfacing(folder,time(end)*60*fps-binSliding/2*60*fps:time(end)*60*fps) = nan;
    
    
    %% DISTRIBUTION FUNC. SINGLE DIVES ------------------------------------- %%
   
    tdive{folder}.range = 0:tbins:240;
    tdive{folder}.hist = histc(allDive,tdive{folder}.range);
    tdive{folder}.pdf = histc(allDive,tdive{folder}.range)/length(allDive)/tbins;
    tdive{folder}.cdf = cumsum(tdive{folder}.pdf*tbins);
    
    % distribution of dive intervals
    intdive{folder}.range = 0:tbins*20:1800;
    intdive{folder}.hist = histc(allIntDive,intdive{folder}.range);
    intdive{folder}.pdf = histc(allIntDive,intdive{folder}.range)/length(allIntDive)/(tbins*20);
    intdive{folder}.cdf = cumsum(intdive{folder}.pdf*tbins*20);
    
    AllSingleDiveT = cat(1,AllSingleDiveT,allDive');
    
%     tdig{folder}.range = 0:tbins:600;
%     tdig{folder}.hist = hist(allDig,tdig{folder}.range);
%     
%     tsurf{folder}.range = 0:tbins:600;
%     tsurf{folder}.hist = hist(allSurf,tsurf{folder}.range);

  %% Distribution of max. dive depth ------------------------------------- %
    DiveDepth{folder}.range = 0:zbins:12;
    DiveDepth{folder}.hist = histc(maxDiveDepth,DiveDepth{folder}.range);
    DiveDepth{folder}.pdf = histc(maxDiveDepth,DiveDepth{folder}.range)/length(maxDiveDepth)/zbins;
    DiveDepth{folder}.cdf = cumsum(DiveDepth{folder}.pdf*zbins);
    
    AllMaxDiveDepth = cat(1,AllMaxDiveDepth,maxDiveDepth');
    
  %% Getting all mean dive speed and work
  meanVelDiveAll = cat(1,meanVelDiveAll,mean(meanVelDive));
  semVelDiveAll = cat(1,semVelDiveAll,std(meanVelDive)/sqrt(length(meanVelDive)));
  workDiveAll = cat(1,workDiveAll,workDive');

    %% Getting all mean surf speed 
  meanVelSurfAll = cat(1,meanVelSurfAll,mean(meanVelSurf));
  semVelSurfAll = cat(1,semVelSurfAll,std(meanVelSurf)/sqrt(length(meanVelSurf)));
    
    %% GET ALL RISK-TAKING INDEX ---
    meanRI = mean(RI);
    semRI = std(RI)/sqrt(length(RI));
    AllRI = cat(1,AllRI,RI');  % concanation of all R
    
    %% Calculate prob. of transition sequence : diving-digging-surfacing
%       ndive = [];
%       cumndive = [];
%       proDiveDigSurf = [];
%       ndive = AllnDive(find(groupLarvae==folder));
%       cumndive = cumsum(ndive);
%       value = logical(allSurfAfDive);
%       for j=1:length(ndive)
%           if j==1
%           proDiveDigSurf(j) = sum(value(1:ndive(j)))/ndive(j);
%           else                    
%           proDiveDigSurf(j) = sum(value(cumndive(j-1)+1:cumndive(j)-1))/ndive(j);
%           end
%       end
%     
        
      
%    overlabIndex = find(ismember(allDive,allBeforeDiveInterval(1,:)));
%    allDigAfDive = allDigAfDive(overlabIndex);
%    allSurfAfDive = allSurfAfDive(overlabIndex);  
      
      
    %% WRITE ALL DATA TO SAVE ONE OUTPUT DATA FILE --------------------- %%
    DnDdata{folder}.Genotype = geno;
    DnDdata{folder}.Condition = cond{folder};
    DnDdata{folder}.TotDiveTime = (totTdive/60/totTime*100);
    DnDdata{folder}.TotDigTime = (totTdig/60/totTime*100);
    DnDdata{folder}.TotSurfTime = (totTsurf/60/totTime*100);
    DnDdata{folder}.DiveProbEachLarva = DiveProb;
    DnDdata{folder}.MeanDiveProb = MeanDiveProb(folder,:);
    DnDdata{folder}.SemDiveProb = SemDiveProb(folder,:);
    DnDdata{folder}.DiveProbAll = pDiving(folder,:);
    DnDdata{folder}.DigProbAll = pDigging(folder,:);
    DnDdata{folder}.SurfProbAll = pSurfacing(folder,:);
    DnDdata{folder}.DiveFreq = (nDive/max(time));
    DnDdata{folder}.MeanDiveFreq = nanmean(nDive/max(time));
    DnDdata{folder}.SemDiveFreq = nanstd(nDive/max(time))/sqrt(nLarvae(folder));
    DnDdata{folder}.MaxDiveDepth = maxDiveDepth;
    DnDdata{folder}.SingleDiveTime = allDive;
    DnDdata{folder}.DiveInterval = allIntDive;
    DnDdata{folder}.DiveTimeDiveInterval = allBeforeDiveInterval; % 1st row: singleDiveTime, 2nd row: intervalBeteenDives
    DnDdata{folder}.IniTimeDive = allIniDive;
    DnDdata{folder}.CDFsingleDiveTimeRange = tdive{folder}.range;
    DnDdata{folder}.CDFsingleDiveTime = tdive{folder}.cdf;
    DnDdata{folder}.CDFmaxDiveDepthRange = DiveDepth{folder}.range;
    DnDdata{folder}.CDFmaxDiveDepth = DiveDepth{folder}.cdf;
    DnDdata{folder}.MeanVelDive = meanVelDive;
    DnDdata{folder}.MeanVelSurf = meanVelSurf;
    DnDdata{folder}.WorkDive = workDive;
    DnDdata{folder}.RiskTakingIndex = RI;
    DnDdata{folder}.PreferenceIndex = PI;
    DnDdata{folder}.MeanRiskTakingIndex = meanRI;
    DnDdata{folder}.SemRiskTakingIndex = semRI;
    DnDdata{folder}.FirstSurfTimeAfterDive = allSurfAfDive;
    DnDdata{folder}.FirstDigTimeAfterDive = allDigAfDive;
%    DnDdata{folder}.ProbTransDiveDigSurf = proDiveDigSurf;
    
    %% plot transition prob. matrix
%     fig2 = figure(2);
%     clim = [0,0.5];
%     subplot(2,1,folder)
%     imagesc(trMatrix{folder},clim)
%     colorbar
%     axis equal
%     title(cond{folder})

%     subplot(2,2,3)
%     imagesc(trMatrix{2},clim)
%     colorbar
%     axis equal
    
    %% PLOT MEAN DIVE PROB. OVER TIME ------------------------------------- %%
    fig6=figure(6);
%     plot(DiveRange,DnDdata{folder}.MeanDiveProb,strcat('o',cc{folder}))
%     hold on
%     plot(DiveRange,DnDdata{folder}.MeanDiveProb+DnDdata{folder}.SemDiveProb,strcat('-',cc{folder}),'LineWidth',1)
%     plot(DiveRange,DnDdata{folder}.MeanDiveProb-DnDdata{folder}.SemDiveProb,strcat('-',cc{folder}),'LineWidth',1)
%     xlabel('Time (min)');
%     ylabel('Probability');
    subplot(2,2,1)
    title('Probability distribution of diving')
    %plot(DiveRange+binSliding/2,DnDdata{folder}.DiveProbAll,strcat('-',cc{folder}),'LineWidth',2)
    plot(DiveRange,pDiving(folder,:),strcat('-',cc{folder}),'LineWidth',2)
    hold on
    %xlim([0 30])
    ylim([0 1])
    xlabel('Time (min)','FontSize',12);
    ylabel('Probability','FontSize',12);
    legend(Labels, 'Location','NorthEast');
    
    subplot(2,2,2)
    title('Probability distribution of digging')
    plot(DiveRange,pDigging(folder,:),strcat('-',cc{folder}),'LineWidth',2)
    hold on
    ylim([0 1])
    xlabel('Time (min)','FontSize',12);
    ylabel('Probability','FontSize',12);
    legend(Labels, 'Location','NorthEast');
    
    subplot(2,2,3)
    title('Probability distribution of surfacing')
    plot(DiveRange,pSurfacing(folder,:),strcat('-',cc{folder}),'LineWidth',2)
    hold on
    ylim([0 1])
    xlabel('Time (min)','FontSize',12);
    ylabel('Probability','FontSize',12);
    legend(Labels, 'Location','NorthEast');
    
    subplot(2,2,4)
    title('Preference index')
    xRange = [];
    xRange = folder + randn(size(PI))*0.1;
    medianPI = nanmedian(PI);
    plot(xRange,PI,'ok','MarkerSize',4)
    xlim([0 maxFolder+1])
    ylim([-1.5 1.5])
    hold on
    line([folder-0.2 folder+0.2],[medianPI medianPI],'color','r','Linewidth',2)
    if (folder==maxFolder)
    set(gca,'xTick',[1:1:maxFolder])
    set(gca, 'XTickLabel', Labels, 'FontSize',8)
    %rotateXLabels(gca(),30);
    end

    
    %% SCATTER PLOTS OF DIGGING AND SURFACING TIME AFTER DIVE ---------- %%
%     fig8 = figure(8);
%     
%     subplot(2,3,3*folder-2)
%     plot(allBeforeDiveInterval(1,:)/60,allDigAfDive/60,strcat('o','b'),'MarkerSize',5)
%     %lsline
%     hold on
%     %plot(nanmean(allDive),nanmean(allDigAfDive),'ok','MarkerSize',5,'MarkerFaceColor','k')
%     xlim([-0.2 240/60])
%     ylim([-0.2 480/60])
%     axis equal
%     title('Digging time after a dive')
%     xlabel('Single dive time (min)');
%     ylabel('Digging time after a dive (min)');
%     
%     
%     
%     subplot(2,3,3*folder-1)
%     %plot(allDive(allSurfAfDive~=0),allSurfAfDive(allSurfAfDive~=0),strcat('o','r'),'MarkerSize',5)
%     plot(allBeforeDiveInterval(1,:)/60,allSurfAfDive/60,strcat('o','r'),'MarkerSize',5)
%     %semilogx(allDive,allSurfAfDive,strcat('o','r'),'MarkerSize',5)
%     %lsline
%     %plot(allMdive,allIntDive,strcat('o',cc{folder}),'MarkerSize',5)
% %    plot(allMdive,allIntDive,'o')
%     hold on
%     %plot(nanmean(allDive),nanmean(allSurfAfDive),'ok','MarkerSize',5,'MarkerFaceColor','k')
%     xlim([-0.2 240/60])
%     ylim([-0.2 480/60])
%     axis equal
%     title('Surfacing time after a dive')
%     xlabel('Single dive time (min)');
%     ylabel('Surfacing time after a dive (min)');    
%     
%     
%         subplot(2,3,3*folder)
%     %plot(allDive(allSurfAfDive~=0),allSurfAfDive(allSurfAfDive~=0),strcat('o','r'),'MarkerSize',5)
%     plot(allBeforeDiveInterval(1,:)/60,allBeforeDiveInterval(2,:)/60,strcat('o','k'),'MarkerSize',5)
%     %semilogx(allDive,allSurfAfDive,strcat('o','r'),'MarkerSize',5)
%     %lsline
%     %plot(allMdive,allIntDive,strcat('o',cc{folder}),'MarkerSize',5)
% %    plot(allMdive,allIntDive,'o')
%     hold on
%     %plot(nanmean(allDive),nanmean(allSurfAfDive),'ok','MarkerSize',5,'MarkerFaceColor','k')
%     xlim([-0.2 240/60])
%     ylim([-0.2 480/60])
%     axis equal
%     title('Time interval bwt. dives')
%     xlabel('Single dive time (min)');
%     ylabel('Time interval between dives (min)');
% 
% [rr,pp] = corr(allBeforeDiveInterval(1,:)',allDigAfDive','type','spearman')
% [rr,pp] = corr(allBeforeDiveInterval(1,:)',allSurfAfDive','type','spearman')
% [rr,pp] = corr(allBeforeDiveInterval(1,:)',allBeforeDiveInterval(2,:)','type','spearman')

    
    %% PLOTING DIVE PARAMETERS ----------------------------------------- %%
    fig5=figure(5);
    subplot(2,2,1)
    plot(DiveDepth{folder}.range, DiveDepth{folder}.cdf, strcat('-',cc{folder}),'LineWidth',2);
    hold on
    legend(Labels, 'Location','SouthEast');
    xlim([3 12])
    ylim([0 1.2])
    xlabel('Max. dive depth (mm)');
    ylabel('Probability');
    
    subplot(2,2,2)
    plot(tdive{folder}.range, tdive{folder}.cdf, strcat('-',cc{folder}),'LineWidth',2);
    hold on
    legend(Labels, 'Location','SouthEast');
    xlim([0 240])
    ylim([0 1.2])
    %plot(digz{folder}.range, digz{folder}.kspdf,strcat('red'));
    xlabel('Single dive time (s)');
    ylabel('Probability');
    
    subplot(2,2,3)
    plot(allDive,maxDiveDepth,strcat('o',cc{folder}),'MarkerSize',5)
    hold on
    legend(Labels, 'Location','SouthEast');
    lsline
    xlim([0 220])
    ylim([0 14])
    xlabel('Sigle dive time (s)');
    ylabel('Max. dive depth (mm)');
    
    subplot(2,2,4)
    plot(intdive{folder}.range, intdive{folder}.cdf, strcat('-',cc{folder}),'LineWidth',2);
    hold on
    legend(Labels, 'Location','SouthEast');
    xlim([0 1800])
    ylim([0 1.2])
    %plot(digz{folder}.range, digz{folder}.kspdf,strcat('red'));
    xlabel('Time interval between dives (s)');
    ylabel('Probability');
    
    
end

% plot Dive frequency   ----------------------------%
% subplot(2,2,4)
% bar(meanFdive,0.5)
% colormap([0.4 0.4 1]);
% ylabel('Dive frequency (1/min)');
% set(gca,'xtick',1:maxFolder)
% set(gca, 'XTickLabel', Labels)
% hold on
% errorbar(1:maxFolder,meanFdive, semFdive,'k','linestyle','none')
% -------

cd(DirHome)

if (EthogramON == 1)
%%  PLOT ETHOGRAM ------------------------------------------------------ %%
ang=0:0.01:2*pi; 
xCir = 20;
yCir = [10 60 110];
lArrow = 8;
gap = 4;

% cCode{1} = 'green';
% cCode{2} = 'blue';
% cCode{3} = 'red';
scal = 0.0;
cCode(1,:) = [scal 1 scal];
cCode(2,:) = [scal 0 1-scal];
cCode(3,:) = [1 scal scal];

% sumTR = zeros(3,3);
% for i=1:maxFolder
%     sumTR = sumTR + trMatrix{i};
% end

fig2 = figure(2);

for i = 1:maxFolder
    subplot(2,6,i)
    title(cond{i})
    xlim([0 40])
    ylim([0 120])
    
    hold on
    axis equal
    axis off
    for k=1:3
        for l=1:3
            wArrow = 1000*trMatrix{i}(k,l); % width of arrow
            cArrow = mean([yCir(k),yCir(l)]);
            yArrow = cArrow-lArrow/2+lArrow*[0 1 1 1+wArrow/lArrow*0.3*sqrt(3) 1 1 0];               % Y-coords of arrow edge points
            xArrow = xCir+wArrow*[-0.1 -0.1 -0.3 0 0.3 0.1 0.1];  % x-coords of arrow edge points
            
            if (k-l==-1)
                xArrow = xArrow - gap;
                hArrow = fill(xArrow,yArrow,'k');    % Plot a red arrow
                text(xCir-22,cArrow,num2str(trMatrix{i}(k,l),'%.3f'),'FontSize',8,'FontWeight', 'bold')
            elseif (k-l==1)
                xArrow = xArrow + gap;  % x-coords of arrow edge points
                [xArrow,yArrow] = rotateData(xArrow,yArrow,xCir+gap,cArrow,180);
                hArrow = fill(xArrow,yArrow,'k');    % Plot a red arrow
                text(xCir+7,cArrow,num2str(trMatrix{i}(k,l),'%.3f'),'FontSize',8,'FontWeight', 'bold')
            elseif (k-l==0)
                %radi = 20*trMatrix{i}(k,k)/(trace(sumTR));
                radi = 20*trMatrix{i}(k,k)/sum(sum(trMatrix{i}));
                xp = radi*cos(ang);
                yp = radi*sin(ang);      
                fill(xCir+xp,yCir(k)+yp,cCode(k,:),'FaceAlpha',0.6)
                text(xCir-7,cArrow-(radi+5)*(k~=2),num2str(trMatrix{i}(k,l),'%.3f'),'FontSize',8,'FontWeight', 'bold')
            end
        end
    end   
end

end

%% PLOT FREQUENCY OF DIVE ---------------------------------------------- %%
fig3=figure(3);
subplot(2,2,1)
bar(meanFdive,0.5)
colormap([0.4 0.4 1]);
ylabel('Dive frequency (1/min)');
ylim([0 0.5])
%set(gca,'xtick',1:maxFolder)
set(gca, 'XTickLabel', Labels, 'FontSize',8)
hold on
errorbar(1:maxFolder,meanFdive, semFdive,'k','linestyle','none')
%rotateXLabels(gca(),30);

subplot(2,2,2)
bar(meanIntDive,0.5)
colormap([0.4 0.4 1]);
ylabel('Time interval btw. dives (s)');
ylim([0 500])
%set(gca,'xtick',1:maxFolder)
set(gca, 'XTickLabel', Labels, 'FontSize',8)
hold on
errorbar(1:maxFolder,meanIntDive,semIntDive,'k','linestyle','none')
%rotateXLabels(gca(),30);

subplot(2,2,3)
bar(meanVelDiveAll,0.5)
colormap([0.4 0.4 1]);
ylabel('Mean Dive Speed (mm/s)');
ylim([0 0.5])
%set(gca,'xtick',1:maxFolder)
set(gca, 'XTickLabel', Labels, 'FontSize',8)
hold on
errorbar(1:maxFolder,meanVelDiveAll,semVelDiveAll,'k','linestyle','none')
%rotateXLabels(gca(),30);

subplot(2,2,4)
bar(meanIntDive,0.5)
colormap([0.4 0.4 1]);
ylabel('Time interval btw. dives (s)');
ylim([0 500])
%set(gca,'xtick',1:maxFolder)
set(gca, 'XTickLabel', Labels, 'FontSize',8)
hold on
errorbar(1:maxFolder,meanIntDive,semIntDive,'k','linestyle','none')
%rotateXLabels(gca(),30);

% subplot(2,2,4)
% boxplot(workDiveAll, groupDive, 'notch', 'on')
% ylabel(gca,'Work (arb. unit)','FontSize',12)
% ylim([-1 20])
% %set(gca,'xtick',1:maxFolder)
% %set(gca, 'XTickLabel', Labels)
% text_h = findobj(gca, 'Type', 'text');
%     
%     for cnt = 1:length(text_h)
%         set(text_h(cnt),    'FontSize', fontSize,...
%                             'Rotation', rotation, ...
%                             'String', Labels{length(Labels)-cnt+1}, ...
%                             'HorizontalAlignment', 'right')
%     end

%% BOXPLOTING % TIME IN EACH BEHAVIORAL MODE --------------------------- %%
if (BoxplotON ==1 )  
    
fig4=figure(4);
subplot(2,3,1)
%boxplot(AlltotTdive,groupLarvae,'notch','on')
boxplot(AlltotTdive,groupLarvae,'notch','off')
title(quan{1})
%ylim([-2.5 42.5])
ylim([-2.5 42.5])
%ylim([-2.5 52.5])
%ylim([-5 105])
set(gca,'YTick',[0:20:100])
ylabel(gca,'% time in behavioural mode','FontSize',12)
%xlabel('Substrate hardness (% agarose)')
% set(gca,'xtick',1:maxFolder)
% set(gca,'XTickLabel',Labels)


text_h = findobj(gca, 'Type', 'text');
    
    for cnt = 1:length(text_h)
        set(text_h(cnt),    'FontSize', fontSize,...
                            'Rotation', rotation, ...
                            'String', Labels{length(Labels)-cnt+1}, ...
                            'HorizontalAlignment', 'right')
    end


subplot(2,3,2)
boxplot(AlltotTdig,groupLarvae,'notch','off')
title(quan{2})
ylim([-5 105])
set(gca,'YTick',[0:20:100])
ylabel(gca,'% time in behavioural mode','FontSize',12)
%xlabel('Substrate hardness (% agarose)')
% set(gca,'xtick',1:maxFolder)
% set(gca,'XTickLabel',Labels)

text_h = findobj(gca, 'Type', 'text');
    
    for cnt = 1:length(text_h)
        set(text_h(cnt),    'FontSize', fontSize,...
                            'Rotation', rotation, ...
                            'String', Labels{length(Labels)-cnt+1}, ...
                            'HorizontalAlignment', 'right')
    end

subplot(2,3,3)
boxplot(AlltotTsurf,groupLarvae,'notch','off')
title(quan{3})
ylim([-5 105])
%ylim([-2.5 52.5])
%ylim([-2.5 42.5])
set(gca,'YTick',[0:20:100])
ylabel(gca,'% time in behavioural mode','FontSize',12)
%xlabel('Substrate hardness (% agarose)')
% set(gca,'xtick',1:maxFolder)
% set(gca,'XTickLabel',Labels)

text_h = findobj(gca, 'Type', 'text');
    
    for cnt = 1:length(text_h)
        set(text_h(cnt),    'FontSize', fontSize,...
                            'Rotation', rotation, ...
                            'String', Labels{length(Labels)-cnt+1}, ...
                            'HorizontalAlignment', 'right')
    end

subplot(2,3,4)
boxplot(AllRI,groupLarvae,'notch','on')
title(quan{4})
%ylim([-0.525 0.025])
ylim([-0.025 1.225])
%set(gca,'YTick',[0:20:100])
ylabel(gca,'Risk-taking index','FontSize',12)
%xlabel('Substrate hardness (% agarose)')
% set(gca,'xtick',1:maxFolder)
% set(gca,'XTickLabel',Labels)

text_h = findobj(gca, 'Type', 'text');
    
    for cnt = 1:length(text_h)
        set(text_h(cnt),    'FontSize', fontSize,...
                            'Rotation', rotation, ...
                            'String', Labels{length(Labels)-cnt+1}, ...
                            'HorizontalAlignment', 'right')
    end

% subplot(2,3,5)
% boxplot(AllSingleDiveT,groupDive,'notch','on')
% title('Single Dive Time')
% ylim([-10 300])
% %set(gca,'YTick',[0:10:50])
% ylabel(gca,'Time (s)','FontSize',12)
% %xlabel('Substrate hardness (% agarose)')
% % set(gca,'xtick',1:maxFolder)
% % set(gca,'XTickLabel',Labels)
% 
% text_h = findobj(gca, 'Type', 'text');
%     
%     for cnt = 1:length(text_h)
%         set(text_h(cnt),    'FontSize', fontSize,...
%                             'Rotation', rotation, ...
%                             'String', Labels{length(Labels)-cnt+1}, ...
%                             'HorizontalAlignment', 'right')
%     end
% 
% subplot(2,3,6)
% boxplot(AllMaxDiveDepth,groupDive,'notch','on')
% title('Maximum Dive Depth')
% ylim([0 12])
% %set(gca,'YTick',[0:10:50])
% ylabel(gca,'Max. Dive Depth (mm)','FontSize',12)
% %xlabel('Substrate hardness (% agarose)')
% % set(gca,'xtick',1:maxFolder)
% % set(gca,'XTickLabel',Labels)
% 
% text_h = findobj(gca, 'Type', 'text');
%     
%     for cnt = 1:length(text_h)
%         set(text_h(cnt),    'FontSize', fontSize,...
%                             'Rotation', rotation, ...
%                             'String', Labels{length(Labels)-cnt+1}, ...
%                             'HorizontalAlignment', 'right')
%     end

fig7=figure(7);   
subplot(2,2,1)
%title('Single Dive Time')
hold on
for folder = 1:maxFolder 
xRange = [];
xRange = folder + randn(size(DnDdata{folder}.SingleDiveTime))*0.1;
medianVal = median(DnDdata{folder}.SingleDiveTime/60);
plot(xRange,DnDdata{folder}.SingleDiveTime/60,'ok','MarkerSize',4)
line([folder-0.2 folder+0.2],[medianVal medianVal],'color','r','Linewidth',2)
end

xlim([0 maxFolder+1])
ylim([-0.5 6])
ylabel('Single dive time (min)')
%set(gca, 'xtick',[1:maxFolder])
%set(gca, 'XTickLabel', Labels, 'FontSize',8)    

% xticklabel_rotate([1:maxFolder],45,Labels,'interpreter','none')


subplot(2,2,2)
%title('Maximum Dive Depth')
hold on
for folder = 1:maxFolder 
xRange = [];
xRange = folder + randn(size(DnDdata{folder}.MaxDiveDepth))*0.1;
medianVal = median(DnDdata{folder}.MaxDiveDepth);
plot(xRange,DnDdata{folder}.MaxDiveDepth,'ok','MarkerSize',3)
line([folder-0.2 folder+0.2],[medianVal medianVal],'color','r','Linewidth',2)
end

xlim([0 maxFolder+1])
ylim([2 12])
ylabel('Max. dive depth (mm)')
%set(gca, 'xtick',[1:maxFolder])
%set(gca, 'XTickLabel', Labels, 'FontSize',8)  
    
% xticklabel_rotate([1:maxFolder],45,Labels,'interpreter','none')

end


if (SaveDataON == 1)
    saveas(fig3,strcat(DirExport,'DiveFreqSpeedWork'),'png');
saveas(fig5,strcat(DirExport,'CDFdive'),'png');
saveas(fig7,strcat(DirExport,'singleDive'),'png');

savefile = strcat(geno,'_DnDdata','.mat');
save(strcat(DirExport,savefile), 'DnDdata');

savefile = strcat(geno,'_1',quan{1},'.mat');
AlltotTimeMode = AlltotTdive;
save(strcat(DirExport,savefile), 'AlltotTimeMode','groupLarvae');
savefile = strcat(geno,'_2',quan{2},'.mat');
AlltotTimeMode = AlltotTdig;
save(strcat(DirExport,savefile), 'AlltotTimeMode','groupLarvae');
savefile = strcat(geno,'_3',quan{3},'.mat');
AlltotTimeMode = AlltotTsurf;
save(strcat(DirExport,savefile), 'AlltotTimeMode','groupLarvae');

savefile = strcat(geno,'_DiveFreq','.mat');
save(strcat(DirExport,savefile), 'AllfDive','groupLarvae');

end



if (SaveStatON == 1)
  
%% Statistics and exporting results
% Multiple comparision of dive speed%
for i=1:length(conditions)
    diveSpeedAll = cat(1,diveSpeedAll,DnDdata{i}.MeanVelDive');
end
[pAnova,table,stats] = anova1(diveSpeedAll,groupDive,'off');
%[c,m,h,nms] = multcompare(stats,'alpha',0.01,'ctype','hsd');
pAnova;
cPvalue = 0.05;
combi = nchoosek(1:maxFolder,2);
nPairs = size(combi,1);
    
% exporting statistics %%
f = fopen(strcat(DirExport,'DiveFreqSpeedWork_statistics.txt'),'wt');
%fprintf(f, 'Filename                                 Mean    STD    NRT(T) FRT(1)  FRT(2)  FRT(3)  SD(FRT)\n'); 

fprintf(f, 'Information of analysis : %s\n',date);
% fprintf(f, 'Agarose: 0.4 percent \n');
% fprintf(f, 'Assaytime: 30 min \n');
% fprintf(f, 'Pre-diffusion time for odor: 20 min \n');

% comparing dive speed
    for i=1:nchoosek(maxFolder,2)
        [h(i),p(i)]=ttest2(DnDdata{combi(i,1)}.MeanVelDive,DnDdata{combi(i,2)}.MeanVelDive,cPvalue/nPairs);
    end
pCorrect = p*nPairs; % bonferroni correction

for i=1:length(conditions)
    fprintf(f, '%s(%i,n= %i)\n',cond{i},i,stats.n(i));
end
fprintf(f, '\n\n');

fprintf(f, 'Multiple comparison of dive speed\n');
fprintf(f, 'One-way ANOVA : P = %1.10e\n',pAnova);
fprintf(f, 'Post hoc : bonferroni@ P < %f\n\n',cPvalue);

for j=1:nchoosek(groupDive(end),2)
    fprintf(f, 'Comparison of dive speed: %s vs. %s\n', cond{combi(j,1)},cond{combi(j,2)});
    fprintf(f, '\n');
    fprintf(f, 'Corrected P-value: %1.4e, Diffenent?: %i,    Uncorrected P-value: %1.4e\n', pCorrect(j), h(j), p(j));
    fprintf(f, '\n\n');
end

% comparing works
    for i=1:nchoosek(maxFolder,2)
        [p(i),h(i)]=ranksum(workDiveAll(groupDive==combi(i,1)),workDiveAll(groupDive==combi(i,2)),cPvalue/nPairs);
    end
pCorrect = p*nPairs;     % bonferroni correction
for j=1:nchoosek(groupDive(end),2)
    fprintf(f, 'Comparison of dive work: %s vs. %s\n', cond{combi(j,1)},cond{combi(j,2)});
    fprintf(f, '\n');
    fprintf(f, 'Corrected P-value: %1.4e, Diffenent?: %i,    Uncorrected P-value: %1.4e\n', pCorrect(j), h(j), p(j));
    fprintf(f, '\n\n');
end

% comparing dive frequency
    for i=1:nchoosek(maxFolder,2)
        [h(i),p(i)]=ttest2(DnDdata{combi(i,1)}.DiveFreq,DnDdata{combi(i,2)}.DiveFreq,cPvalue/nPairs);
    end
pCorrect = p*nPairs; % bonferroni correction

for j=1:nchoosek(maxFolder,2)
    fprintf(f, 'Comparison of dive frequency: %s vs. %s\n', cond{combi(j,1)},cond{combi(j,2)});
    fprintf(f, '\n');
    fprintf(f, 'Corrected P-value: %1.4e, Diffenent?: %i,    Uncorrected P-value: %1.4e\n', pCorrect(j), h(j), p(j));
    fprintf(f, '\n\n');
end

% comparing dive time intervals
    for i=1:nchoosek(maxFolder,2)
        [h(i),p(i)]=ttest2(DnDdata{combi(i,1)}.DiveInterval,DnDdata{combi(i,2)}.DiveInterval,cPvalue/nPairs);
    end
pCorrect = p*nPairs; % bonferroni correction

for j=1:nchoosek(groupDive(end),2)
    fprintf(f, 'Comparison of dive time intervals: %s vs. %s\n', cond{combi(j,1)},cond{combi(j,2)});
    fprintf(f, '\n');
    fprintf(f, 'Corrected P-value: %1.4e, Diffenent?: %i,    Uncorrected P-value: %1.4e\n', pCorrect(j), h(j), p(j));
    fprintf(f, '\n\n');
end

fclose(f); 

% compare % time in behaviral mode
% ranksum test
for j=1:3
    
    if (j==1) %disp('TotalDivingTime')
    for i=1:nchoosek(maxFolder,2)
        [p(j,i),h(j,i)]=ranksum(AlltotTdive(groupLarvae==combi(i,1)),AlltotTdive(groupLarvae==combi(i,2)));
    end
    end
    
    if (j==2) %disp('TotalDiggingTime')
    for i=1:nchoosek(maxFolder,2)
        [p(j,i),h(j,i)]=ranksum(AlltotTdig(groupLarvae==combi(i,1)),AlltotTdig(groupLarvae==combi(i,2)));
    end
    end
    
    if (j==3) %disp('TotalSurfacingTime')
    for i=1:nchoosek(maxFolder,2)
        [p(j,i),h(j,i)]=ranksum(AlltotTsurf(groupLarvae==combi(i,1)),AlltotTsurf(groupLarvae==combi(i,2)));
    end
    end
           
end
% exporting results
f = fopen(strcat(DirExport,'statistics_TimeSpentModes.txt'),'wt');
%fprintf(f, 'Filename                                 Mean    STD    NRT(T) FRT(1)  FRT(2)  FRT(3)  SD(FRT)\n'); 

fprintf(f, 'Information of analysis : %s\n',date);
fprintf(f, 'Agarose: 0.4 percent \n');
fprintf(f, 'Assaytime: %i min \n',analTime);
fprintf(f, 'thDive: %1.4f \n',thDive);
%fprintf(f, 'Pre-diffusion time for odor: 20 min \n');

for i=1:maxFolder
    fprintf(f, '%s(%i): %i\n',cond{i}, i, nLarvae(i));
end
fprintf(f, '\n\n');

for i=1:nchoosek(maxFolder,2)
    fprintf(f, 'Comparison of total time in behavioral mode: %s vs. %s\n', cond{combi(i,1)},cond{combi(i,2)});
    fprintf(f, '\n');
    for j=1:3
        fprintf(f, '%s: %1.8f %i\n', quan{j}, p(j,i), h(j,i));
    end
    fprintf(f, '\n');
end

% comparing risk-taking index
    for i=1:nchoosek(maxFolder,2)
        [p(i),h(i)]=ranksum(AllRI(groupLarvae==combi(i,1)),AllRI(groupLarvae==combi(i,2)),cPvalue/nPairs);
    end
pCorrect = p*nPairs;     % bonferroni correction
for j=1:nchoosek(groupDive(end),2)
    fprintf(f, 'Comparison of risk-taking index: %s vs. %s\n', cond{combi(j,1)},cond{combi(j,2)});
    fprintf(f, '\n');
    fprintf(f, 'Corrected P-value: %1.4e, Diffenent?: %i,    Uncorrected P-value: %1.4e\n', pCorrect(j), h(j), p(j));
    fprintf(f, '\n\n');
end

% comparing single dive times
    for i=1:nchoosek(maxFolder,2)
        [p(i),h(i)]=ranksum(AllSingleDiveT(groupDive==combi(i,1)),AllSingleDiveT(groupDive==combi(i,2)),cPvalue/nPairs);
    end
pCorrect = p*nPairs;     % bonferroni correction
for j=1:nchoosek(groupDive(end),2)
    fprintf(f, 'Comparison of single dive times: %s vs. %s\n', cond{combi(j,1)},cond{combi(j,2)});
    fprintf(f, '\n');
    fprintf(f, 'Corrected P-value: %1.4e, Diffenent?: %i,    Uncorrected P-value: %1.4e\n', pCorrect(j), h(j), p(j));
    fprintf(f, '\n\n');
end

% comparing max. dive depths
    for i=1:nchoosek(maxFolder,2)
        [p(i),h(i)]=ranksum(AllMaxDiveDepth(groupDive==combi(i,1)),AllMaxDiveDepth(groupDive==combi(i,2)),cPvalue/nPairs);
    end
pCorrect = p*nPairs;     % bonferroni correction
for j=1:nchoosek(groupDive(end),2)
    fprintf(f, 'Comparison of max. dive depths: %s vs. %s\n', cond{combi(j,1)},cond{combi(j,2)});
    fprintf(f, '\n');
    fprintf(f, 'Corrected P-value: %1.4e, Diffenent?: %i,    Uncorrected P-value: %1.4e\n', pCorrect(j), h(j), p(j));
    fprintf(f, '\n\n');
end

fclose(f);


%% SAVE RESULTS AND FIGURES -------------------------------------------- %%
if (RasterON == 1)
    saveas(fig1,strcat(DirExport,'Rasterplot'),'png');
end
if (EthogramON == 1)
    saveas(fig2,strcat(DirExport,'Ethogram'),'png');
end
if (BoxplotON == 1)
    saveas(fig4,strcat(DirExport,'Boxplot'),'png');
end

end

figure
bar(meanVelSurfAll,0.5)
ylabel('Mean surfacing speed (mm/s)');
ylim([0 0.5])
%set(gca,'xtick',1:maxFolder)
set(gca, 'XTickLabel', Labels, 'FontSize',8)
hold on
errorbar(1:maxFolder,meanVelSurfAll,semVelSurfAll,'k','linestyle','none')
% rotateXLabels(gca(),30);

