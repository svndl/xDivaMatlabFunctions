function [screenImage, preludeImage] = generateStimset(timingXDiva, video, stimset)
       
    nBarPairs = 4;
    nBars = 2*nBarPairs;
    maxDisp = max(stimset.fig.dispSteps/video.pix2arcmin);
    switch stimset.stimGeom
        case 'Hbars'
            limits.x = [-maxDisp*0.5 + stimset.dotSizePix, maxDisp*0.5 + video.width_pix - stimset.dotSizePix];
            limits.y = [-maxDisp*0.5 + stimset.dotSizePix, maxDisp*0.5 + video.height_pix/nBars - stimset.dotSizePix];
            dx = 0;
            dy = video.height_pix/nBars;
            shift.x = 0;
            shift.y = .5*video.height_pix;
        case 'Vbars'
            limits.x = [-maxDisp*0.5 + stimset.dotSizePix, maxDisp*0.5 + video.width_pix/nBars - stimset.dotSizePix];
            limits.y = [-maxDisp*0.5 + stimset.dotSizePix, maxDisp*0.5 + video.height_pix - stimset.dotSizePix];
            dx = video.width_pix/nBars;
            dy = 0;
            shift.x = .5*video.width_pix;
            shift.y = 0; 
    end    
    % start position for all dots all sweep steps
    switch stimset.viewMode
        case 'Fullscreen'
            % don't care
        case 'Square'
            minDim = min(video.height_pix, video.width_pix);
            video.height_pix = minDim;
            video.width_pix = minDim;
    end
    barInfo.modEnv = stimset.modEnv;
    barInfo.modType = stimset.modType;
    barInfo.numDots = round(stimset.numDots/nBars);    
    
    
    timing.nSteps = timingXDiva.nCoreSteps;
    timing.nFrames = getUniqueFrames(stimset.modEnv, timingXDiva);    
    timing.framesPerStep = timing.nFrames/timing.nSteps;
    %% bgr bar
    bgrInfo = barInfo;
    
    bgrInfo.dispSteps = stimset.fig2bgr*stimset.fig.dispSteps/video.pix2arcmin;
    bgrInfo.cohSteps = stimset.bgr.cohSteps*ones(1, timing.nSteps)/100;
    bgrInfo.corrSteps = stimset.bgr.corrSteps*ones(1, timing.nSteps);
    
    
    %% figure bar    
    figInfo = barInfo;
    
    figInfo.dispSteps = stimset.fig.dispSteps/video.pix2arcmin;
    figInfo.cohSteps = stimset.fig.cohSteps/100;
    figInfo.corrSteps = stimset.fig.corrSteps;
    
       
    % width of one bar is dy (or dx)
    dots.L.x = [];
    dots.L.y = [];
    dots.R.x = [];
    dots.R.y = [];
    
    LimitsBgr = limits;
    
    for b = 1:nBarPairs
     
        limitsFig.x = LimitsBgr.x + dx;
        limitsFig.y = LimitsBgr.y + dy;

        startPosBgr.L = getStartPos(barInfo.numDots, 1, LimitsBgr);
        startPosBgr.R = startPosBgr.L;
        
        startPosFig.L = getStartPos(barInfo.numDots, 1, limitsFig);
        startPosFig.R = startPosFig.L;
        
        dotsBgr = mkBarFrames(bgrInfo, startPosBgr, LimitsBgr, timing);    
        dotsFig = mkBarFrames(figInfo, startPosFig, limitsFig, timing);
        dots.L.x = cat(1, dots.L.x, dotsBgr.L.x, dotsFig.L.x);    
        dots.L.y = cat(1, dots.L.y, dotsBgr.L.y, dotsFig.L.y);    
        dots.R.x = cat(1, dots.R.x, dotsBgr.R.x, dotsFig.R.x);    
        dots.R.y = cat(1, dots.R.y, dotsBgr.R.y, dotsFig.R.y);    
        
        LimitsBgr.x = LimitsBgr.x + 2*dx;
        LimitsBgr.y = LimitsBgr.y + 2*dy;
    end
    
    %% draw fixation
    stimset.drawCross = 0;
    if (strcmp(stimset.fixPoint, 'Cross'))
        stimset.drawCross = 1;
    end
    
    screenImage = DrawDotFrames(video, stimset, dots);
        
    %% generate prelude
    [preludeInfo, timingPrelude] = generatePrelude(timingXDiva, stimset);
    preludeImage = [];
    if (~isempty(preludeInfo))
        limitsPrelude.x = [0 video.width_pix];
        limitsPrelude.y = [0 video.height_pix];
    
        startPosPrelude.L = getStartPos(preludeInfo.numDots, timingPrelude.nFrames, limitsPrelude);
        startPosPrelude.R = startPosPrelude.L;
        preludeInfo.drawCross = 0;
        preludeImage = DrawDotFrames(video, preludeInfo, startPosPrelude);
    end
    
    %% throw away G (2nd dim)
    screenImage(:, :, 2, :) = [];
    preludeImage(:, :, 2, :) = [];
end
%% generate one bar
function dots = mkBarFrames(barInfo, startPos, limits, timing)
    dots.R.x = [];
    dots.R.y = [];
    dots.L.x = [];
    dots.L.y = [];
    dotColors = [];   
        
    [lrSign, dotShift] = mkShift(timing.nFrames, barInfo.dispSteps, ...
        barInfo.modType, barInfo.modEnv);
    framesPerStep = timing.nFrames/timing.nSteps;
    for ns = 1:timing.nSteps
        % generate dot frames (observing coherence and ) 
        ch = barInfo.cohSteps(ns);
        cr = barInfo.corrSteps(ns);
       
        [dotFrames, colorsSweep] = mkStepSweepFrames(startPos, framesPerStep, limits, ch, cr);
        
        dotShift_s.x = dotShift.x(:, (ns - 1)*framesPerStep + 1:ns*framesPerStep);
        dotShift_s.y = dotShift.y(:, (ns - 1)*framesPerStep + 1:ns*framesPerStep);

        shifted = shiftDots(dotFrames, dotShift_s, lrSign, limits, 1);

        
        dots.L.x = [dots.L.x shifted.L.x];
        dots.L.y = [dots.L.y shifted.L.y];
       
        dots.R.x = [dots.R.x shifted.R.x];
        dots.R.y = [dots.R.y shifted.R.y];
       
%         startPos.L.x = shifted.L.x(:, end);
%         startPos.L.y = shifted.L.y(:, end);
%        
%         startPos.R.x = shifted.R.x(:, end);
%         startPos.R.y = shifted.R.y(:, end);
       
        %dotColors = [dotColors repmat(colorsSweep,  [1 timing.figFramesPerCycle])];
    end
end
%% draw dots
function screenImage = DrawDotFrames(video, stimset, dots)
    %% OPEN PTB
    whichScreen = max(Screen('Screens'));
    myRect = [0 0 video.height_pix video.width_pix];
    PsychImaging('PrepareConfiguration');
    Screen('Preference', 'SkipSyncTests', 2);
    Screen('Preference', 'Verbosity', 0);
    Screen('Preference', 'VisualDebugLevel', 0);

   
    [window, windowRect] = PsychImaging('OpenWindow', whichScreen, 0, myRect, 24, [], [], video.aaFactor); % 8-bit/channel precision
    Screen(window, 'BlendFunction', GL_ONE, GL_ONE); % allows overlapping dot regions to blend properly, GL_ONE = PTB default var

    %% test flips
   
    for ntestFlips = 1:20
        vbl = Screen('Flip', window);
        ifi = Screen('GetFlipInterval', window);
    end
    idxUpdate = 1;
    waitframes = 0.5;
    nDotsUpdate = size(dots.R.x, 2);
    screenImage = zeros(video.height_pix, video.width_pix, 3, nDotsUpdate, 'int8');    
    tic
    
    %% fix cross settings
    cross.dimPix = 30;
    cross.x = [-cross.dimPix cross.dimPix 0 0];
    cross.y = [0 0 -cross.dimPix cross.dimPix];
    cross.width_pix = 2;
    cross.color = .8*WhiteIndex(whichScreen);
    [cross.x0, cross.y0] = RectCenter(windowRect);  
    
    while(idxUpdate <= nDotsUpdate)
        lDots = [dots.L.x(:, idxUpdate), dots.L.y(:, idxUpdate)]';
        rDots = [dots.R.x(:, idxUpdate), dots.R.y(:, idxUpdate)]';
        
        
        %scrLCenter = [video.width_pix/2 video.height_pix/2];
        scrLCenter = [0 0];       
        scrRCenter = scrLCenter;

        Screen('DrawDots', window, lDots, stimset.dotSizePix, [255 0 0], scrLCenter, 2); % 2, blended circles
        Screen('DrawDots', window, rDots, stimset.dotSizePix, [0 0 255], scrRCenter, 2); % 2, blended circles
        
        if (stimset.drawCross)
            Screen('DrawLines', window, [cross.x; cross.y],...
                cross.width_pix, cross.color, [cross.x0, cross.y0], 2);
        end
        
        vbl = Screen('Flip', window, vbl + waitframes*ifi);
       
        screenImage(:, :, :, idxUpdate) = Screen('GetImage', window );
        idxUpdate = idxUpdate + 1;
    end
    toc
    Screen('DrawingFinished', window, 0);
    Screen('Flip', window);
    Screen('CloseAll');
    sca;
end
%%
function out = getStartPos(nDots, nFrames, limits)
%     out.x = limits.x(2)*(2*rand(nDots, nFrames) - 1) + limits.x(1);
%     out.y = limits.y(2)*(2*rand(nDots, nFrames) - 1) + limits.y(1);
    out.x = (limits.x(2) - limits.x(1))*rand(nDots, nFrames) + limits.x(1);
    out.y = (limits.y(2) - limits.y(1))*rand(nDots, nFrames) + limits.y(1);
end

function [dotFrames, dotColorsSweep] = mkStepSweepFrames(startPos, nFrames, limits, ch, cr)
    % total number of dots
    nDots = size(startPos.L.x, 1);

    nCohDots = round(nDots*ch);
    nuCohDots = nDots - nCohDots;
    
    nCohCorrDots = round(nCohDots*abs(cr));
    nCohuCorrDots = nCohDots - nCohCorrDots;
    
    nuCohCorrDots = round(nuCohDots*abs(cr));    
    nuCohuCorrDots = nDots - nCohDots - nuCohCorrDots;
    
    % generate four types of dot motion
    
    %% coherent correlated dots
    CohCorr.L.x = repmat(startPos.L.x(1:nCohCorrDots), [1 nFrames]);
    CohCorr.L.y = repmat(startPos.L.y(1:nCohCorrDots), [1 nFrames]);
    
    CohCorr.R.x = repmat(startPos.L.x(1:nCohCorrDots), [1 nFrames]);
    CohCorr.R.y = repmat(startPos.L.y(1:nCohCorrDots), [1 nFrames]);
    
    %% coherent uncorrelated dots
    % left eye gets the start mat 
    CohnCorr.L.x = repmat(startPos.L.x(nCohCorrDots + 1:nCohCorrDots + nCohuCorrDots), [1 nFrames]);
    CohnCorr.L.y = repmat(startPos.L.y(nCohCorrDots + 1:nCohCorrDots + nCohuCorrDots), [1 nFrames]);
    
    % right eye will get different dots
    rightStartPos = getStartPos(nCohuCorrDots, 1, limits);
    CohnCorr.R.x = repmat(rightStartPos.x, [1 nFrames]);
    CohnCorr.R.y = repmat(rightStartPos.y, [1 nFrames]);
    
    %% incoherent correlated dots
    IncohStartPos = getStartPos(nuCohCorrDots, nFrames - 1, limits);
    
    nCohCorr.L.x = [startPos.L.x(nCohCorrDots + 1 + nCohuCorrDots : ...
        nCohCorrDots + nCohuCorrDots + nuCohCorrDots)  IncohStartPos.x];
    nCohCorr.L.y = [startPos.L.y(nCohCorrDots + 1 + nCohuCorrDots : ...
        nCohCorrDots + nCohuCorrDots + nuCohCorrDots) IncohStartPos.y];
    
    nCohCorr.R.x = [startPos.L.x(nCohCorrDots + 1 + nCohuCorrDots : ...
        nCohCorrDots + nCohuCorrDots + nuCohCorrDots)  IncohStartPos.x];
    nCohCorr.R.y = [startPos.L.y(nCohCorrDots + 1 + nCohuCorrDots : ...
        nCohCorrDots + nCohuCorrDots + nuCohCorrDots) IncohStartPos.y];

    %% incoherent uncorrelated dots
    
    nCohnCorr.L = getStartPos(nuCohuCorrDots, nFrames, limits);
    nCohnCorr.R = getStartPos(nuCohuCorrDots, nFrames, limits);
    

    dotFrames.L.x = cat(1, CohCorr.L.x, CohnCorr.L.x, nCohCorr.L.x, nCohnCorr.L.x);
    dotFrames.L.y = cat(1, CohCorr.L.y, CohnCorr.L.y, nCohCorr.L.y, nCohnCorr.L.y);
    dotFrames.R.x = cat(1, CohCorr.R.x, CohnCorr.R.x, nCohCorr.R.x, nCohnCorr.R.x);
    dotFrames.R.y = cat(1, CohCorr.R.y, CohnCorr.R.y, nCohCorr.R.y, nCohnCorr.R.y);
    
    %color-code each dot
    dotColorsSweep = cat(1, ones(nCohCorrDots, 1), sign(cr)*ones(nCohuCorrDots, 1), ...
        ones(nuCohCorrDots, 1), sign(cr)*ones(nuCohuCorrDots, 1));
end
%% perform dot shift
function shifted = shiftDots(dotFrames, shift, lrSign, limits, dotFramesPerCycle)
    
    dotShift.x = repmat(shift.x, [size(dotFrames.L.x, 1) 1]);
    dotShift.y = repmat(shift.y, [size(dotFrames.L.x, 1) 1]);
    
    % add shift to coherent dots only
    posLx = limits.x(1) + rem(dotFrames.L.x + dotShift.x, limits.x(2) - limits.x(1));
    posLy = limits.y(1) + rem(dotFrames.L.y + dotShift.y, limits.y(2) - limits.y(1));
    posRx = limits.x(1) + rem(dotFrames.R.x + lrSign*dotShift.x, limits.x(2) - limits.x(1));
    posRy = limits.y(1) + rem(dotFrames.R.y + lrSign*dotShift.y, limits.y(2) - limits.y(1));
    
    
    shifted.L.x = my_repelem( posLx, 1, dotFramesPerCycle);
    shifted.L.y = my_repelem( posLy, 1, dotFramesPerCycle);
    
    shifted.R.x = my_repelem( posRx, 1, dotFramesPerCycle);
    shifted.R.y = my_repelem( posRy, 1, dotFramesPerCycle);
end
%% generate prelute

function [preludeInfo, preludeTiming] = generatePrelude(timingXDiva, stimset)
    if (timingXDiva.nPreludeBins)
       preludeInfo = [];
       preludeTiming = [];
    end
    
    switch timingXDiva.preludeType
        case 'Blank'
            nFrames = 1;
            numDots = 0;
        case 'Dynamic'
            nFrames = timingXDiva.figFramesPerCycle/timingXDiva.dotFramesPerCycle;
            numDots = stimset.numDots;
        case 'Static'
            nFrames = 1;
            numDots = stimset.numDots;
    end
 
    preludeInfo.numDots = numDots;
    preludeInfo.dotSizePix = stimset.dotSizePix;
    preludeInfo.corrSteps = 1;
    preludeInfo.cohSteps = 0;
    preludeInfo.dispSteps = 0;
    
    preludeTiming.nFrames = nFrames;
    preludeTiming.nSteps = timingXDiva.nPreludeBins;
end
%% Number of unique frames depending on the mod env
function nFrames = getUniqueFrames(modEnv, timingXDiva)
    nFrames = timingXDiva.nCoreSteps*timingXDiva.figFramesPerCycle/...
        timingXDiva.dotFramesPerCycle;

    switch modEnv
        case 'motion on-off'
            nFrames = nFrames*2;
    end
end
