clear all;close all;
delete(gcp);
numcells=8;

numsweeps=10;
%%%%%%%%%%%parameters for detectPSPs%%%%%%%%%%%%%%%%%%%%



readParams.startTime = 1;
readParams.endTime = 9.7;
% readParams.endTime = 4.7;
readParams.highCut = 5000;
readParams.roundTimesVal = 1;
readParams.numThreshSweeps = 10;
readParams.doDownSample = 1;

readParams.doPlotSweeps = true;
readParams.doMarkEvents = false;
readParams.doPlotEvents = false;

detectParams.timesRMS = 2;

... minimum allowable amplitude for alpha functions (in units of samples)
    detectParams.minAmp = -1000;
... maximum allowable amplitude
    detectParams.maxAmp = -5;
... minimum allowable tau for alpha functions (in units of samples)
    detectParams.minTau = 0;%original 500e-6,50e-6
... maximum allowable tau
    detectParams.maxTau = inf;%original 5e-3,worked 50e-3
... minimum allowable yOffset for alpha functions (in units of mV)
    detectParams.minYOffset = -2000;
... maximum allowable yOffset
    detectParams.maxYOffset = 2000;
... minimum allowable decay tau
    detectParams.minDecay = 0; % CHANGE FROM 1E-3 ON 11-05-19 (Scott had set to 3e-3)
... maximum allowable decay tau
    detectParams.maxDecay = inf; % CHANGE FROM 50E-3 ON 11-05-19  (Scott had set to 20e-3)
... threshold used to determine if the change of derivative is great
    ... enough to separate out two alpha functions
    detectParams.derThresh = 5;%original 10,5 worked for cell4, derthresh=5 for 3 and 4
... second EPSP is not recognized if it is closer than this value to
    ... the previous one (in units of samples)
    detectParams.closestEPSPs = 1E-3;
% threshold for standard error above which a fit with multiple alphas is
% attempted
detectParams.errThresh = 0.08;
% 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
detectParams.dataFilterType = 3;
% 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
detectParams.derFilterType = 3;
% length of data filter
detectParams.dataFilterLength = 5E-3;
% length of derivative filter
detectParams.derFilterLength = 2E-3;
% if set to 1 then debugging figures appear
detectParams.debugging = 0;
% index of first data point
detectParams.dataStart = 1;
% forces a graphical output even if other outputs are taken
detectParams.forceDisplay = 0;
% turns off best fitting of start time and decay Tau when 1
detectParams.noFit = 0;

detectParams.threshVal=-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



i=16;
results={};

filepath='/home/ben/Desktop/5ht3a_project/abf_files/';
%set up a parallel pool so we can use speedy parallel loops; here we use 4
%workers
parpool('local',4);
tic,
for j=1:3
    disp(strcat('condition:',num2str(j)));
    ppp=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The following if/else statements are particular to my files, you
    %should modify them so they correclty call your files%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==1
        detectParams.derThresh=10;
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/3_18_15/15318';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/3_18_15/cell1_31815_results.mat';
        if j==1
            range=0:9;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=10:19;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=20:29;
            obslength(i,j)=length(range)*10;
        end
        
    elseif i==2
        detectParams.derThresh=8;
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/3_19_15/15319';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/3_19_15/cell1_31915_results.mat';
        if j==1
            range=10:19;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=20:29;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=30:39;
            obslength(i,j)=length(range)*10;
        end
        
    elseif i==3
        detectParams.derThresh=8;
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/3_19_15/15319';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/3_19_15/cell1_31915_10uM_results.mat';
        if j==1
            range=40:49;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=50:59;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=60:69;
            obslength(i,j)=length(range)*10;
        end
        
    elseif i==4
        detectParams.derThresh=8;
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/3_19_15/15319';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/3_19_15/cell1_31915_100uM_results.mat';
        if j==1
            range=70:79;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=80:89;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=90:99;
            obslength(i,j)=length(range)*10;
        end
        
    elseif i==5%between 1 and 10 uM
        detectParams.derThresh=15;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/3_20_15/15320';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/3_20_15/cell1_32015_results_thresh15.mat';
        if j==1
            range=80:89;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=90:99;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=100:109;
            obslength(i,j)=length(range)*10;
        end
        
    elseif i==6%1 uM
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/3_20_15/15320';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/3_20_15/cell1_32015_cell2_results.mat';
        if j==1
            range=0:9;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=10:19;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=20:29;
            obslength(i,j)=length(range)*10;
            
        end
        
    elseif i==7%10 uM
        detectParams.derThresh=8;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/3_20_15/15320';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/3_20_15/32015_cell3_results.mat';
        if j==1
            range=30:39;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=40:49;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=50:59;
            obslength(i,j)=length(range)*10;
            
        end
        
    elseif i==8%1 uM
        
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/3_16_15/15316';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/3_16_15/31615_cell3_results.mat';
        if j==1
            range=0:9;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=10:19;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=20:29;
            obslength(i,j)=length(range)*10;
            
        end
        
    elseif i==9
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/3_20_15/15320';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/3_20_15/cell1_32015_100uM_results.mat';
        if j==1
            range=50:59;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=60:69;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=70:79;
            obslength(i,j)=length(range)*10;
            
        end
        
    elseif i==10
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_3_15/15403';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_3_15/cell1_40315_1uM_results.mat';
        if j==1
            range=0:9;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=10:19;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=20:29;
            obslength(i,j)=length(range)*10;
            
        end
        
    elseif i==11
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_3_15/15403';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_3_15/cell2_40315_10uM_results.mat';
        if j==1
            range=51:55;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=56:65;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=66:70;
            obslength(i,j)=length(range)*10;
            
        end
    elseif i==13
        detectParams.derThresh=15;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_13_15/15413';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_13_15/cell1_41315_1uM_results.mat';
        if j==1
            range=0:9;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=10:19;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=20:29;
            obslength(i,j)=length(range)*10;
            
        end
    elseif i==14
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_17_15/15417';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_17_15/cell2_41715_100uM_results.mat';
        if j==1
            range=67:76;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=77:86;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=87:96;
            obslength(i,j)=length(range)*10;
            
        end
    elseif i==15
        setthresh=1;
        detectParams.derThresh=15;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_17_15/15417';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_17_15/cell1_41715_100uM_results.mat';
        if j==1
            range=37:46;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=47:56;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=57:66;
            obslength(i,j)=length(range)*10;
            
        end
    elseif i==16
        setthresh=0;
        detectParams.derThresh=5;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_21_15/15421';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_21_15/cell1_42115_100uM_resultsderthresh5.mat';
        if j==1
            range=9:18;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=19:28;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=29:38;
            obslength(i,j)=length(range)*10;
            
        end
        
    elseif i==17
        setthresh=0;
        detectParams.derThresh=5;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_23_15/15423';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_23_15/cell1_42315_100uM_resultsderthresh5.mat';
        if j==1
            range=0:9;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=10:19;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=20:35;
            obslength(i,j)=length(range)*10;
            
        end
        
    elseif i==18
        setthresh=0;
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_23_15/15423';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_23_15/cell2_42315_100uM_results.mat';
        if j==1
            range=50:59;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=60:69;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=70:85;
            obslength(i,j)=length(range)*10;
            
        end
                
    elseif i==19
        setthresh=0;
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_24_15/15424';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_24_15/cell1_42415_100uM_results.mat';
        if j==1
            range=0:9;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=10:19;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=20:35;
            obslength(i,j)=length(range)*10;
            
        end
                        
    elseif i==20
        setthresh=0;
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_24_15/15424';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_24_15/cell2_42415_100uM_results.mat';
        if j==1
            range=38:47;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=48:57;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=58:72;
            obslength(i,j)=length(range)*10;
            
        end
                                
    elseif i==21
        setthresh=0;
        detectParams.derThresh=10;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_24_15/15424';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_24_15/cell3_42415_100uM_results.mat';
        if j==1
            range=73:82;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=83:92;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=93:107;
            obslength(i,j)=length(range)*10;
            
        end
        
    elseif i==22
        setthresh=1;
        detectParams.derThresh=5;%original 10
        PSPsDown=1;
        filepath2='/home/ben/Desktop/5ht3a_project/abf_files/4_27_15/15427';
        ressave='/home/ben/Desktop/5ht3a_project/abf_files/4_27_15/cell1_42715_100uM_resultsderthresh5.mat';
        if j==1
            range=0:9;
            obslength(i,j)=length(range)*10;
        elseif j==2
            range=10:19;
            obslength(i,j)=length(range)*10;
            
        elseif j==3
            range=20:35;
            obslength(i,j)=length(range)*10;
            
        end
    end
    
    for k=1:length(range)
        
        if range(k)<10
            filepath=strcat(filepath2,'00',num2str(range(k)),'.abf');
        elseif range(k)>9&&range(k)<100
            filepath=strcat(filepath2,'0',num2str(range(k)),'.abf');
        elseif range(k)>99
            filepath=strcat(filepath2,num2str(range(k)),'.abf');
        end
        
        x=abfload(filepath);
        x=x(1.2e3:end,1,:);
        
        resultstemp={};
        %parallel for loop for speed; note how to set up a parallel pool
        %properly with parpool earlier in the code
        parfor kk=1:numsweepsl
            
            disp(strcat('sweep:',num2str((k-1)*10+kk)));
            SF=10000;
            downSampleFreq=SF;
            %filter data with a butterworth filter with a low pass cutoff below 1000
            %Hz
            [b,a]=butter(6,1000/SF,'low');
            %xx=filtfilt(b,a,x(8e3:end,1,:));
            
            
            
            filtData=filtfilt(b,a,x(:,1,kk));
            %     results=cat(1,results,detectPSPs(xx(:,1,i)));
            %the following if/else statements either call detectPSPs and
            %force it to plot (setthresh=1) so you can see if your
            %detection criteria are adequate or run through all the data and save it (which
            %should only be done after you've properly set your detection
            %params)
            if setthresh
                detectPSPs(filtData, PSPsDown, ...
                    'minAmp', detectParams.minAmp, ...
                    ... minimum allowable amplitude for alpha functions (in units of samples)
                    'maxAmp', detectParams.threshVal, ...
                    ... maximum allowable amplitude
                    'minTau', (detectParams.minTau * downSampleFreq), ...
                    ... minimum allowable tau for alpha functions (in units of samples)
                    'maxTau', (detectParams.maxTau * downSampleFreq), ...
                    ... maximum allowable tau
                    'minYOffset', detectParams.minYOffset, ...
                    ... minimum allowable yOffset for alpha functions (in units of mV)
                    'maxYOffset', detectParams.maxYOffset, ...
                    ... maximum allowable yOffset
                    'minDecay', (detectParams.minDecay * downSampleFreq), ...
                    ... minimum allowable decay tau
                    'maxDecay', (detectParams.maxDecay * downSampleFreq), ...
                    ... maximum allowable decay tau
                    'derThresh', (detectParams.derThresh), ...
                    ... threshold used to determine if the change of derivative is ...
                    ... great enough to separate out two alpha functions
                    'closestEPSPs', (detectParams.closestEPSPs * downSampleFreq), ...
                    ... second EPSP is not recognized if it is closer than this ...
                    ... value to the previous one (in units of samples)
                    'errThresh', (detectParams.errThresh), ...
                    ... threshold for standard error above which a fit with ...
                    ... multiple alphas is attempted
                    'dataFilterType', detectParams.dataFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'derFilterType', detectParams.derFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'dataFilterLength', (detectParams.dataFilterLength * downSampleFreq), ...
                    ... length of data filter
                    'derFilterLength', (detectParams.derFilterLength * downSampleFreq), ...
                    ... length of derivative filter
                    'debugging', detectParams.debugging, ...
                    ... if set to 1 then debugging figures appear
                    'dataStart', detectParams.dataStart, ...
                    ... index of first data point
                    'forceDisplay', detectParams.forceDisplay,...
                    ... forces a graphical output even if other outputs are taken
                    'noFit', detectParams.noFit)
                pause();
                close;
            else
                resultstemp{kk}=detectPSPs(filtData, PSPsDown, ...
                    'minAmp', detectParams.minAmp, ...
                    ... minimum allowable amplitude for alpha functions (in units of samples)
                    'maxAmp', detectParams.threshVal, ...
                    ... maximum allowable amplitude
                    'minTau', (detectParams.minTau * downSampleFreq), ...
                    ... minimum allowable tau for alpha functions (in units of samples)
                    'maxTau', (detectParams.maxTau * downSampleFreq), ...
                    ... maximum allowable tau
                    'minYOffset', detectParams.minYOffset, ...
                    ... minimum allowable yOffset for alpha functions (in units of mV)
                    'maxYOffset', detectParams.maxYOffset, ...
                    ... maximum allowable yOffset
                    'minDecay', (detectParams.minDecay * downSampleFreq), ...
                    ... minimum allowable decay tau
                    'maxDecay', (detectParams.maxDecay * downSampleFreq), ...
                    ... maximum allowable decay tau
                    'derThresh', (detectParams.derThresh), ...
                    ... threshold used to determine if the change of derivative is ...
                    ... great enough to separate out two alpha functions
                    'closestEPSPs', (detectParams.closestEPSPs * downSampleFreq), ...
                    ... second EPSP is not recognized if it is closer than this ...
                    ... value to the previous one (in units of samples)
                    'errThresh', (detectParams.errThresh), ...
                    ... threshold for standard error above which a fit with ...
                    ... multiple alphas is attempted
                    'dataFilterType', detectParams.dataFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'derFilterType', detectParams.derFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'dataFilterLength', (detectParams.dataFilterLength * downSampleFreq), ...
                    ... length of data filter
                    'derFilterLength', (detectParams.derFilterLength * downSampleFreq), ...
                    ... length of derivative filter
                    'debugging', detectParams.debugging, ...
                    ... if set to 1 then debugging figures appear
                    'dataStart', detectParams.dataStart, ...
                    ... index of first data point
                    'forceDisplay', detectParams.forceDisplay,...
                    ... forces a graphical output even if other outputs are taken
                    'noFit', detectParams.noFit);
                
                %                         pause();
                %                         close;
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            %             ppp=ppp+1;
            ... turns off best fitting of start time and decay Tau when 1
                %                     pause();
            %                 close;
            
        end
        for mmm=1:numsweeps
            results{j,(k-1)*10+mmm}=resultstemp{mmm};
        end
        clear resultstemp;
        
        
        
    end
    
    
    
    
end
toc
delete(gcp);
save(ressave,'results');
%     
    
