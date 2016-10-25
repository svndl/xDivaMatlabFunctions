function [lrSign, shift] = mkShift(nFrames, displacement, modType, modEnv)
    % figure out shift sign for each eye
    shift.x = zeros(1, nFrames);
    shift.y = zeros(1, nFrames);
    nCycles = numel(displacement);
    startMod = zeros(nCycles, 1);
    endMod = zeros(nCycles, 1);
    switch modType
        %% zero->A modulation
        case 'NearZero_x'
            % towards: left dots move right            
            lrSign = -1;
            coord = 'x';
            startMod = displacement;
        case 'FarZero_x'
            % away: left dots move left            
            lrSign = -1;
            coord = 'x';
            startMod = -displacement;
        case 'NearZero_y'
            % towards: left dots move right            
            lrSign = -1;
            coord = 'y';
            startMod = displacement;
        case 'FarZero_y'
            % away: left dots move left            
            lrSign = -1;
            coord = 'y';
            startMod = -displacement;
        case 'LeftZero'
            lrSign = 1;
            coord = 'x';
            startMod = -displacement;
        case 'RightZero'
            lrSign = 1;
            coord = 'x';
            startMod = displacement;
        case 'UpZero'
            lrSign = 1;
            coord = 'y';
            startMod = -displacement;
        case 'DownZero'
            lrSign = 1;
            startMod = displacement;
            coord = 'y';
        %% -A->A modulation
        case 'DownUp'
            lrSign = 1;
            coord = 'y';
            startMod = displacement;
            endMod = -displacement;
        case 'NearFar_x'
            lrSign = -1;
            coord = 'x';
            startMod = displacement;
            endMod = -displacement;
        case 'NearFar_y'
            lrSign = -1;
            coord = 'y';
            startMod = displacement;
            endMod = -displacement;
        case 'RightLeft'
            lrSign = 1;
            coord = 'x';
            startMod = displacement;
            endMod = -displacement;            
    end
    
    %% generating the shift values for the whole sweep
    
    nFramesPerCycle = nFrames/nCycles;
    cycle1 = round(.5*nFramesPerCycle);
    cycle2 = nFramesPerCycle - cycle1;
    displacement = [];    
    switch modEnv
        case 'square'
            for d = 1:nCycles
            	displacement = [displacement repmat(startMod(d), [1 cycle1]) repmat(endMod(d), [1 cycle2])];
            end
        case 'triangular'
            shiftVal = .5*(startMod(1) + endMod(1));            
            for d = 1:nCycles
                d = linspace
                displacement = [];
            end
        case 'sinusoidal'
            shiftD = displacement*sin(linspace(-pi/2*(startD/displacement), ...
                pi/2*(endD/displacement), nFrames));
        case 'sawtooth'
            for d = 1:nCycles           
            end
        case 'motion on-off'
            LastPos = endMod(1);
            onoffStart = endMod;
            onoffEnd = startMod;
            quarterCycle = round(.25*nFramesPerCycle);

            displacement = [];
            for d = 1:nCycles
                q1 = linspace(LastPos, onoffEnd(d), quarterCycle + 1);
                q2 = onoffEnd(d)*ones(1, quarterCycle);
                q3 = linspace(onoffEnd(d), onoffStart(d), quarterCycle + 1);
                q4 = onoffStart(d)*ones(1, quarterCycle);
                displacement = [displacement q1(2:end) q2 q3(2:end) q4];
                LastPos = q4(end);
            end
    end
    shift.(coord) = displacement;    
end 