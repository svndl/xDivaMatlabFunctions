function pmf_RandomDotsStereoMotion( varargin )
%
%   Shows a correlated figure region which potentially modulates in depth
%   over a background whose correlation level can be any value from 0-1.
%
%	xDiva Matlab Function paradigm

%	"pmf_" prefix is not strictly necessary, but helps to identify
%	functions that are intended for this purpose.

%   Each Matlab Function paradigm will have its own version of this file

%	First argument is string for selecting which subfunction to use
%	Additional arguments as needed for each subfunction

%	Each subfunction must conclude with call to "assignin( 'base', 'output', ... )",
%	where value assigned to "output" is a variable or cell array containing variables
%	that xDiva needs to complete the desired task.

    if nargin > 0, aSubFxnName = varargin{1}; else error( 'pmf_RandomDotsStereoMotion.m called with no arguments' ); end

% these next three are shared by nested functions below, so we create
% them in this outermost enclosing scope.
    definitions = MakeDefinitions;
    parameters = {};
    timing = {};
    videoMode = {};

    % some useful functional closures...
    CFOUF = @(varargin) cellfun( varargin{:}, 'uniformoutput', false );
    AFOUF = @(varargin) arrayfun( varargin{:}, 'uniformoutput', false );
    
    %lambda functions
    PVal = @( iPart, x ) ParamValue( num2str(iPart), x );
    PVal_S = @(x) ParamValue( 'S', x );
        
    
    
    % PVal_B = @(x) ParamValue( 'B', x );
    % PVal_1 = @(x) ParamValue( 1, x );
    % PVal_2 = @(x) ParamValue( 2, x );

    aaFactor = 4; 
    maxDispAmin = 70;   
    ppath = setPathParadigm;
    screenRedGunData = load(fullfile(ppath.info, 'SonyTV_RedGunValues.txt'));

    try
        switch aSubFxnName
            case 'GetDefinitions', GetDefinitions;
            case 'ValidateParameters', ValidateParameters;
            case 'MakeMovie', MakeMovie;
        end
    catch tME
        errLog = fopen(fullfile(ppath.log, 'ErrorLog.txt'), 'a+');
        display(tME.message);
        for e = 1: numel(tME.stack)
            fprintf(errLog, ' %s ', tME.stack(e).file);
            fprintf(errLog, ' %s ', tME.stack(e).name);        
            fprintf(errLog, ' %d\n', tME.stack(e).line);        
        end
        fclose(errLog);
        rethrow( tME ); % this will be caught by xDiva for runtime alert message
    end

    function rV = ParamValue( aPartName, aParamName )
        % Get values for part,param name strings; e.g "myViewDist = ParamValue( 'S', 'View Dist (cm)' );"
        tPart = parameters{ ismember( { 'S' 'B' '1' '2' }, {aPartName} ) }; % determine susbscript to get {"Standard" "Base" "Part1" "Part2"} part cell from parameters
        rV = tPart{ ismember( tPart(:,1), { aParamName } ), 2 }; % from this part, find the row corresponding to aParamName, and get value from 2nd column
    end
    function rV = GetParamArray( aPartName, aParamName )
        
        % For the given part and parameter name, return an array of values
        % corresponding to the steps in a sweep.  If the requested param is
        % not swept, the array will contain all the same values.
        
        % tSpatFreqSweepValues = GetParamArray( '1', 'Spat Freq (cpd)' );
        
        % Here's an example of sweep type specs...
        %
        % definitions{end-2} =
        % 	{
        % 		'Fixed'         'constant'   { }
        % 		'Contrast'      'increasing' { { '1' 'Contrast (pct)' } { '2' 'Contrast (pct)' } }
        % 		'Spat Freq'      'increasing' { { '1' 'Spat Freq (cpd)' } { '2' 'Spat Freq (cpd)' } }
        % 	}
        
        T_Val = @(x) timing{ ismember( timing(:,1), {x} ), 2 }; % get the value of timing parameter "x"
        tNCStps = T_Val('nmbCoreSteps');
        tSweepType = PVal_S('Sweep Type');
        
        % we need to construct a swept array if any of the {name,value} in definitions{5}{:,3}
        
        [ ~, tSS ] = ismember( tSweepType, definitions{end-2}(:,1) ); % the row subscript in definitions{5} corresponding to requested sweep type
        % determine if any definitions{5}{ tSS, { {part,param}... } } match arguments tPartName, tParamName
        IsPartAndParamMatch = @(x) all( ismember( { aPartName, aParamName }, x ) );
        tIsSwept = any( cellfun( IsPartAndParamMatch, definitions{end-2}{tSS,3} ) ); % will be false for "'Fixed' 'constant' { }"
        
        if ~tIsSwept
            rV = ones( tNCStps, 1 ) * ParamValue(  aPartName, aParamName );
        else
            tStepType = PVal_S('Step Type');
            tIsStepLin = strcmpi( tStepType, 'Lin Stair' );
            tSweepStart = PVal_S('Sweep Start');
            tSweepEnd = PVal_S('Sweep End');
            if tIsStepLin
                rV = linspace( tSweepStart, tSweepEnd, tNCStps )';
            else
                rV = logspace( log10(tSweepStart), log10(tSweepEnd), tNCStps )';
            end
        end
        
    end

    function rV = MakeDefinitions
        % for "ValidateDefinition"
        % - currently implementing 'integer', 'double', 'nominal'
        % - types of the cells in each parameter row
        % - only "standard" type names can be used in the "type" fields
        % - 'nominal' params should have
        %       (a) at least one item
        %		(b) value within the size of the array
        % - all other params (so far) should have empty arrays of items
        
        rV = { ...
            
            % - Parameters in part_S must use standard parameter names
            % - 'Sweep Type' : at least 'Fixed' sweep type must be defined as first item in list
            % - 'Modulation' : at least 'None' modulation type must be defined
            % - 'Step Type'  : at least 'Lin Stair' type must be defined,
            %                  first 4 step types are reserved, custom step types can only be added after them
        
            % "Standard" part parameters - common to all paradigms, do not modify names.
            {
                'View Dist (cm)'	100.0	'double' {}
                'Mean Lum (cd)'     1.0	'double' {} % under default calibration, this is (0,0,0)
                'Fix Point'         'None'	'nominal' { 'None' 'Cross'}
                'Sweep Type'        'Disp Amp'	'nominal' { 'Fixed', 'Disp Amp', 'Fig Correlation', 'Fig Coherence'}
                'Step Type'         'Lin Stair'	'nominal' { 'Lin Stair' 'Log Stair' }
                'Sweep Start'       0.0	'double' {}
                'Sweep End'         8.0 'double' {}
                'Modulation'        'NearZero_x' 'nominal' { 'None' 'NearZero_x' 'FarZero_x' 'NearFar_x' ...
                'NearZero_y' 'FarZero_y' 'NearFar_y' ...
                'RightZero' 'LeftZero' 'RightLeft'...
                'DownZero' 'UpZero' 'DownUp'}
            }
            % 3D motions
            % Modulation type (left/right eye) 
            % ZeroNear 0->towards
            % ZeroFar 0->away
            % NearFar towards->away
            % 2D motions 
            % ZeroRight -> 0->right
            % ZeroLeft 0->left
            % RightLeft right->left
            % None (same as Questions: coh dots static, incoherent dots depend on lifetime settings (agreement)) 
       
        
            % "Base" part parameters - paradigm specific parameters that apply to unmodulated parts of the stimulus
            {
                'ModInfo'			8.0 			'double'	{} % 'Displacement (amin)'
                'Mod Envelope'      'square' 		'nominal'	{'square' 'motion on-off'}
                'Stereo'			'Red-Blue'		'nominal'	{'Red-Blue' 'Top-Bottom'}
                'Geometry'			'Hbars' 		'nominal'	{'Hbars' 'Vbars'}
                'Stimulus Extent'	'Fullscreen'	'nominal'	{'Fullscreen' 'Square'} % 'SqrOnBlack' 'SqrOnMean' could be included as memory intense as 'Fullscreen'
                'Spatial Frequency' 10.0            'double'	{}
            }
        
            % "Part1" - parameters that apply to part of stimulus that carries first frequency tag.
            % "Cycle Frames" must be first parameter
            {
                'Cycle Frames'	   30.0 'integer'	{}  % framerate(Hz)/stimFreq(Hz) Modulation frequency every other frame (2 Hz) 
                'notused'           0.0	'double'	{}	% 'notused' params will be set to 0 by xDiva
                'notused'           0.0	'double'	{}	% 'notused' params will be set to 0 by xDiva
                'notused'           0.0	'double'	{}	% 'notused' params will be set to 0 by xDiva
                'Fig Corr (-1:1)'	1.0	'double'	{}  % L-R eye correlation (dot position + color: -1 same dots colors opposite; 1 same dots same colors ) for FIGURE
                'Fig Coh (0:100)' 100.0	'double'	{}  % L-R motion coherence (rand/directional) Figure
                'Bgr Corr (-1:1)'	1.0	'double'	{}  % L-R eye correlation (dot position + color: -1 same dots colors opposite; 1 same dots same colors ) for Background
                'Bgr Coh (0:100)' 100.0	'double'	{}  % L-R motion coherence (rand/directional) background
                'Bgr RelAmp (-1:1)'	1.0	'double'	{}  % L-R motion coorelation (direction of motion wrt 'figure' motion) 
            }
            % Ampl = 0 -> static coherent; incoherent (lifetime = 1 frame)
            % => boiling motion; incoherent (inf lifetime) => static.
         
        
            %% <Questions>
            % 1. Lifetime (for coherent dots) = inf (confirm yes/no)
            % 2. Lifetime (for incoherent dots) 
            %       options: 
            %           inf (brownian motion) => static if Amp = 0.
            %           1 frame => boiling if Amp = 0. => SELECTED
        
            %% </Questions>
        
        
            % "Part2" - parameters that apply to part of stimulus that carries second frequency tag.
            {
                'Cycle Frames'			3.0	'integer'	{} % Update frequency (every 3rd frame -> 20 HZ)
                'Contrast (pct)'	  100.0 'double'	{}
                'Diameter (amin)'		1.0 'double'	{}
                'Density (1/deg^2)'		5.0 'double'	{}
            }
        
            % Sweepable parameters
            % The cell array must contain as many rows as there are supported Sweep Types
            % 1st column (Sweep Types) contains Sweep Type as string
            % 2nd column (Stimulus Visiblity) contains one of the following strings,
            % indicating how stimulus visibility changes when corresponding swept parameter value increases:
            %   'constant' - stimulus visibility stays constant
            %   'increasing' - stimulus visibility increases
            %   'decreasing' - stimulus visibility decreases
            % 3rd column contains a single-row cell array of pairs, where each pair is a single-row cell
            % array of 2 strings: { Part name, Parameter name }
        
            % If sweep affects only one part, then you only need one
            % {part,param} pair; if it affects both parts, then you need both
            % pairs, e.g. for "Contrast" and "Spat Freq" below
        
            {
                'Fixed'				'constant'   { }
                'Disp Amp'          'increasing' { { 'B' 'ModInfo' } }
                'Fig Correlation'  	'increasing' { { '1' 'Fig Corr (-1:1)' } }
                'Fig Coherence'		'increasing' { { '1' 'Fig Coh (0:100)' } }
            }
        
            % ModInfo information
            % The cell array must contain as many rows as there are supported Modulations
            % 1st column (Modulation) contains one of the supported Modulation typs as string
            % 2nd column contains the name of the ModInfo parameter as string
            % 3rd column (default value) contains default value of the ModInfo
            % parameter for this Modulation
            {
        
                'None'          'ModInfo'          0.0
                'NearZero_x'    'Disp Amp (amin)' 20.0
                'FarZero_x'     'Disp Amp (amin)' 20.0
                'NearFar_x'     'Disp Amp (amin)' 20.0
                'NearZero_y'	'Disp Amp (amin)' 20.0
                'FarZero_y'     'Disp Amp (amin)' 20.0
                'NearFar_y'     'Disp Amp (amin)' 20.0
                'RightZero'     'Disp Amp (amin)' 20.0
                'LeftZero'      'Disp Amp (amin)' 20.0
                'RightLeft'     'Disp Amp (amin)' 20.0 
                'DownZero'      'Disp Amp (amin)' 20.0
                'UpZero'        'Disp Amp (amin)' 20.0
                'DownUp'        'Disp Amp (amin)' 20.0                
            }
            % check direction    
            % Required by xDiva, but not by Matlab Function
            {
                'Version'					1
                'Adjustable'				true
                'Needs Unique Stimuli'		false % ###HAMILTON for generating new stimuli every time
                'Supports Interleaving'		false
                'Part Name'                 { 'Pattern' 'Dots'}
                'Frame Rate Divisor'		{ 2 1 } % {even # frames/cycle only, allows for odd-- makes sense for dot update}
                'Max Cycle Frames'			{ 120 6 } % i.e. -> 0.5 Hz, 10 Hz
                'Allow Static Part'			{ true true }
            }
        };
    end

    function GetDefinitions
        assignin( 'base', 'output', MakeDefinitions );
    end

    function ValidateParameters
        % xDiva invokes Matlab Engine command:
        
        % pmf_<subParadigmName>( 'ValidateParameters', parameters, timing, videoMode );
        % "parameters" here is an input argument. Its cellarray hass the
        % same structure as "defaultParameters" but each parameter row has only first two
        % elements
        
        % The "timing" and "videoMode" cellarrays have the same row
        % structure with each row having a "name" and "value" elements.
        
        
        [ parameters, timing, videoMode ] = deal( varargin{2:4} );
        
        sweepType = PVal('S','Sweep Type');
        isSwept = ~strcmp(sweepType,'Fixed');
        
        VMVal = @(x) videoMode{ ismember( videoMode(:,1), {x} ), 2 };
        
        width_pix = VMVal('widthPix');
        height_pix = VMVal('heightPix');

        width_cm = VMVal('imageWidthCm');
        viewDistCm = PVal('S','View Dist (cm)');
        width_deg = 2 * atand( (width_cm/2)/viewDistCm );
        pix2arcmin = ( width_deg * 60 ) / width_pix;
        height_cm = VMVal('imageHeightCm');
        height_deg = 2 * atand( (height_cm/2)/viewDistCm );
        pix2arcmin_h = ( height_deg * 60 ) / height_pix;
        
        validationMessages = {};
        
        ValidateDotSize;
        ValidateDisplacementParams;
        ValidateContrast;
        ValidateModulation;
        ValidateModulationEnvelope;
        ValidateSizeSettings;
        ValidateUpdateFrequency;
        ValidateStimFrequency;
        
        parametersValid = isempty( validationMessages );
        
        % _VV_ Note the standard 'output' variable name
        output = { parametersValid, parameters, validationMessages };
        assignin( 'base', 'output', output );
        
        function CorrectParam( aPart, aParam, aVal )
            tPartLSS = ismember( { 'S' 'B' '1' '2' }, {aPart} );
            tParamLSS = ismember( parameters{ tPartLSS }(:,1), {aParam} );
            parameters{ tPartLSS }{ tParamLSS, 2 } = aVal;
        end
        
        function AppendVMs(aStr), validationMessages = cat(1,validationMessages,{aStr}); end
        
        function ValidateStimFrequency
            extent = PVal('B', 'Stimulus Extent');
            geom = PVal('B', 'Geometry');
            size = PVal('B', 'Spatial Frequency');
            switch extent
                case 'Square'
                    amin2pix = min(width_pix, height_pix)/min(60*width_deg, 60*height_deg);
                case 'Fullscreen'
                    switch geom
                        case 'Hbars'
                           amin2pix = 1/pix2arcmin;
                        case 'Vbars'
                            amin2pix = 1/pix2arcmin_h;
                    end
            end
               
            stimSizePix = amin2pix*size;
            if (rem(stimSizePix, 1)>0)
                stimSizeAmin = round(stimSizePix)/amin2pix;
                CorrectParam('B', 'Spatial Frequency', stimSizeAmin);                               
                AppendVMs(sprintf(...
                    'Requested spatial frequencyis not an integer number of pixels, correcting to nearest possible value: %3.4f amin.',...
                    stimSizeAmin));
            end
                       
        end
        function ValidateUpdateFrequency
            figFramesPerCycle = PVal(1,'Cycle Frames'); 
            dotFramesPerCycle = PVal(2,'Cycle Frames');
            if rem(.5*figFramesPerCycle, dotFramesPerCycle)
                AppendVMs( sprintf('Invalid Dot Update Setting, half of cycle frames should be a multiple of dot update frames') );
            end                
        end
        
        function ValidateSizeSettings
            sizeSetting = PVal('B', 'Spatial Frequency');
            stimGeom = PVal('B', 'Geometry');
            
            if ~strcmp(stimGeom,'Hbars') && ~strcmp(stimGeom,'Vbars') && ~strcmp(sizeSetting,'default')
                CorrectParam( 'B', 'Spatial Frequency', 'default' );
                AppendVMs( sprintf('Invalid Size Setting. ''default'' is the only allowed setting for Geometry setting ''%s''',stimGeom) );
            end
        end
        
        function ValidateModulation
            if strcmp(PVal('S','Modulation'),'None')
                CorrectParam('S','Modulation','ZeroNear');
                AppendVMs('Modulation choice ''None'' is not supported. Resetting to an allowed choice, ''ZeroNear''.');
            end
        end
        
       function ValidateModulationEnvelope
            modEnv = PVal('B', 'Mod Envelope');
            
            if strcmp(modEnv, 'motion on-off')
                framesPerStep = timing{ismember(timing(:, 1), 'nmbFramesPerStep'), 2};
                framesPerCycle = PVal(1,'Cycle Frames');
                Fig2StepRatio = framesPerStep/framesPerCycle;
                if (mod(Fig2StepRatio, 2) > 0)
                    %% TODO
                    AppendVMs(sprintf('Modulation envelope ''motion off-on'' must have odd ratio %f', Fig2StepRatio));
                end
           end
        end
       
        function ValidateDotSize
            % For PTB code to work properly, dot size must be a non-zero
            % integer number of pixels
            
            % convert dot diameter in amin into pixels:
            dotSizeAmin = PVal('2','Diameter (amin)');
            dotSizePix = dotSizeAmin/pix2arcmin;
            
            if mod(dotSizePix,1) %decimal
                dotSizePix = round(dotSizePix);
                if dotSizePix < 1
                    dotSizePix = 1;
                end
                dotSizeAmin = dotSizePix*pix2arcmin;
                CorrectParam('2', 'Diameter (amin)', dotSizeAmin);
                AppendVMs(sprintf(...
                    'Invalid Dot diam., which must be an integer number of pixels, corrected to nearest possible value: %3.4f amin.',...
                    dotSizeAmin));
            end
            if dotSizePix < 1
                dotSizePix = 1;
                dotSizeAmin = dotSizePix*pix2arcmin;
                CorrectParam('2', 'Diameter (amin)', dotSizeAmin);
                AppendVMs(sprintf(...
                    'Invalid Dot diam., which must be at least 1 pixel, corrected to nearest possible value: %3.4f amin.',...
                    dotSizeAmin));
            end
        end
        
        function ValidateDisplacementParams
            % for now this function only checks if the displacement is too
            % large or too small
            
            minDisp = pix2arcmin/aaFactor;
            
            if isSwept && strcmp(sweepType,'Displacement')
                sweepStart = PVal('S','Sweep Start');
                if ~validDisplacement(sweepStart)
                    correctDisp(sweepStart,'Sweep Start','S')
                end
                sweepEnd = PVal('S','Sweep End');
                if ~validDisplacement(sweepEnd)
                    correctDisp(sweepEnd,'Sweep End','S')
                end
            else
                figDisp = PVal('B','ModInfo');
                if ~validDisplacement(figDisp)
                    correctDisp(figDisp,'ModInfo','B')
                end
            end
            
            function isValid = validDisplacement(dispIn)
                if dispIn < minDisp || dispIn > maxDispAmin
                    isValid = false;
                else
                    isValid = true;
                end
            end
            
            function correctDisp(valOld,valName,partName)
                if valOld < minDisp
                    CorrectParam(partName,valName,minDisp);
                    AppendVMs(sprintf('Invalid %s value, which must be at least %3.4f amin. Now corrected to this value.',valName,minDisp));
                elseif valOld > maxDispAmin
                    CorrectParam(partName,valName,maxDispAmin);
                    AppendVMs(sprintf('Invalid %s value, which can be at most %3.4f amin. Now corrected to this value.',valName,maxDispAmin));
                end
            end
            
        end
        function ValidateFigCorr
            % background correlation levels must be on the interval [0,1]
            
            % get bgrCorr values
            if isSwept && strcmp(sweepType,'Bgr Corr')
                sweepStart = PVal('S','Sweep Start');
                if ~validCorr(sweepStart)
                    correctCorr(sweepStart,'Sweep Start','S');
                end
                sweepEnd = PVal('S','Sweep End');
                if ~validCorr(sweepEnd)
                    correctCorr(sweepEnd,'Sweep End','S')
                end
            else
                bgrCorr = PVal('B','Bgr Corr (0-1)');
                if ~validCorr(bgrCorr)
                    correctCorr(bgrCorr,'Bgr Corr (0-1)','B')
                end
            end
            
            function isValid = validCorr(corrIn)
                if corrIn<0 || corrIn>1
                    isValid = false;
                else
                    isValid = true;
                end
            end
            
            function correctCorr(valOld,valName,partName)
                if valOld < 0
                    CorrectParam(partName,valName,0);
                    AppendVMs(sprintf('Invalid %s value, which cannot be less than 0. Now corrected to 0.',valName));
                elseif valOld > 1
                    CorrectParam(partName,valName,1);
                    AppendVMs(sprintf('Invalid %s value, which cannot be more than 1. Now corrected to 1.',valName));
                end
            end
            
        end
        
        function ValidateBgrCorr
            % background correlation levels must be on the interval [0,1]
            
            % get bgrCorr values
            if isSwept && strcmp(sweepType,'Bgr Corr')
                sweepStart = PVal('S','Sweep Start');
                if ~validCorr(sweepStart)
                    correctCorr(sweepStart,'Sweep Start','S');
                end
                sweepEnd = PVal('S','Sweep End');
                if ~validCorr(sweepEnd)
                    correctCorr(sweepEnd,'Sweep End','S')
                end
            else
                bgrCorr = PVal('B','Bgr Corr (0-1)');
                if ~validCorr(bgrCorr)
                    correctCorr(bgrCorr,'Bgr Corr (0-1)','B')
                end
            end
            
            function isValid = validCorr(corrIn)
                if corrIn<0 || corrIn>1
                    isValid = false;
                else
                    isValid = true;
                end
            end
            
            function correctCorr(valOld,valName,partName)
                if valOld < 0
                    CorrectParam(partName,valName,0);
                    AppendVMs(sprintf('Invalid %s value, which cannot be less than 0. Now corrected to 0.',valName));
                elseif valOld > 1
                    CorrectParam(partName,valName,1);
                    AppendVMs(sprintf('Invalid %s value, which cannot be more than 1. Now corrected to 1.',valName));
                end
            end
            
        end
        
        function ValidateContrast
            % should check if requested contrast is not available due to
            % number of color bits on the system. ###
            
            if PVal( 2, 'Contrast (pct)' ) > 100
                CorrectParam( '2', 'Contrast (pct)', 100 );
                AppendVMs( 'Invalid Part 2 Contrast (pct): too high, corrected to 100.' );
            end
            
            if PVal( 2, 'Contrast (pct)' ) < 0
                CorrectParam( '2', 'Contrast (pct)', 0 );
                AppendVMs( 'Invalid Part 2 Contrast (pct): too low, corrected to 0.' );
            end
            
        end
        
    end

    function MakeMovie
        % ---- GRAB & SET PARAMETERS ----
        [ parameters, timing, videoMode, trialNumber ] = deal( varargin{2:5} );
        save('pmf_RandomDotsStereoMotion_MakeMovie.mat', 'parameters', 'timing', 'videoMode');
        TRVal = @(x) timing{ ismember( timing(:,1), {x} ), 2 };
        VMVal = @(x) videoMode{ ismember( videoMode(:,1), {x} ), 2 };
        
        needsUnique = definitions{end}{3,2};
        needsImFiles = true;
        preludeType = {'Dynamic', 'Blank', 'Static'};
        % timing/trial control vars
        stimsetTiming.nCoreSteps = TRVal('nmbCoreSteps');
        stimsetTiming.nCoreBins = TRVal('nmbCoreBins');
        stimsetTiming.nPreludeBins = TRVal('nmbPreludeBins');
        stimsetTiming.framesPerStep = TRVal('nmbFramesPerStep');
        stimsetTiming.framesPerBin = TRVal('nmbFramesPerBin');
        stimsetTiming.preludeType = preludeType{1 + TRVal('preludeType')};
        stimsetTiming.isBlankPrelude = stimsetTiming.preludeType == 1;
        stimsetTiming.nCoreFrames = stimsetTiming.framesPerStep * stimsetTiming.nCoreSteps;
        stimsetTiming.nPreludeFrames = stimsetTiming.nPreludeBins * stimsetTiming.framesPerBin;
        stimsetTiming.nTotalFrames = 2 * stimsetTiming.nPreludeFrames + stimsetTiming.nCoreFrames;
        stimsetTiming.figFramesPerCycle = PVal(1,'Cycle Frames'); % part 1 = Figure
        stimsetTiming.dotFramesPerCycle = PVal(2,'Cycle Frames'); % part 2 = Dots
        stimsetTiming.minUpdate = min([stimsetTiming.figFramesPerCycle, stimsetTiming.dotFramesPerCycle]);
        
        % screen vars
        video.width_pix = VMVal('widthPix');
        video.height_pix = VMVal('heightPix');
        video.width_cm = VMVal('imageWidthCm');
        video.height_cm = VMVal('imageHeightCm');
        video.frameRate = VMVal('nominalFrameRateHz');
        video.viewDistCm = PVal('S','View Dist (cm)');
        video.stimContrast = PVal(2, 'Contrast (pct)')/100;
        video.aaFactor = aaFactor;
        %video.stimExtent = PVal(1,'Stimulus Extent');
        
        video.cm2pix = video.width_pix/video.width_cm;
        video.width_deg = 2 * atand( (video.width_cm/2)/video.viewDistCm );
        video.pix2arcmin = ( video.width_deg * 60 ) / video.width_pix;
        video.width_dva = video.width_pix*(video.pix2arcmin/60);
        video.height_dva = video.height_pix*(video.pix2arcmin/60);
        video.redMod = screenRedGunData(screenRedGunData(:,1)==video.width_pix, 2)/255;
        % value of red gun intensity viewed through red lens (3.0186 cd/m2),
        % ~matched to max intensity blue through blue lens (3.03 cd/m2)
        
        % stim vars
        stimset.sweepType = PVal('S','Sweep Type');
        stimset.isSwept = ~strcmp(stimset.sweepType,'Fixed');
        
        
        stimset.modType = PVal('S', 'Modulation');
        stimset.modEnv = PVal('B', 'Mod Envelope');
        stimset.viewMode = PVal('B', 'Stimulus Extent');
        stimset.stimGeom = PVal('B','Geometry');
        stimset.sizeSetting = PVal('B', 'Spatial Frequency');
        stimset.dotSizeAmin = PVal(2, 'Diameter (amin)');
        stimset.dotDensity = PVal(2, 'Density (1/deg^2)');
        % convert dot density to number of dots for PTB:
        stimset.numDots = round(stimset.dotDensity*(video.width_dva * video.height_dva ));
        % convert dot diameter in amin into pixels for PTB:
        stimset.dotSizePix = stimset.dotSizeAmin/video.pix2arcmin;
        stimset.dotSizeAmin = stimset.dotSizePix*video.pix2arcmin;
        
        stimset.maxDispAmin = maxDispAmin;
        stimset.stimSizePix = video.width_pix*(PVal('B', 'Spatial Frequency')/video.width_deg);
        
        stimset.fig.corrSteps = GetParamArray('1','Fig Corr (-1:1)');
        stimset.fig.dispSteps = GetParamArray('B','ModInfo');
        stimset.fig.cohSteps = GetParamArray('1','Fig Coh (0:100)');
       
        stimset.nStims = PVal('B', 'Spatial Frequency');
        stimset.bgr.corrSteps = PVal('1','Bgr Corr (-1:1)');
        stimset.bgr.cohSteps = PVal('1','Bgr Coh (0:100)');
        stimset.fig2bgr = PVal('1', 'Bgr RelAmp (-1:1)');
        
        stimset.fixPoint = PVal('S', 'Fix Point');
        save('xDivaOutput.mat', 'stimset', 'video', 'stimsetTiming');
        
        [rIms, rPrelude] = generateStimset(stimsetTiming, video, stimset);
        
        % stimset timing
        nUniqueFrames = size(rIms, 4);
        UniqueFramesPerStep = nUniqueFrames/stimsetTiming.nCoreSteps;
        rImSeq = [];
        nCyclesPerStep = stimsetTiming.framesPerBin/(UniqueFramesPerStep*stimsetTiming.dotFramesPerCycle);
        
        for s = 1:stimsetTiming.nCoreSteps
            stepFrames = zeros(stimsetTiming.framesPerBin/nCyclesPerStep, 1);
            stepFrames(1:stimsetTiming.dotFramesPerCycle:end) = (s - 1)*UniqueFramesPerStep + 1:s*UniqueFramesPerStep;
            rImSeq = cat(1, rImSeq, repmat(stepFrames, [nCyclesPerStep 1]));
        end
        
        %prelude-postlude timing
        nUniquePreludeFrames = size(rPrelude, 4);
        rImSeqPrelude = [];
        nCyclesPerStepPrelude = stimsetTiming.framesPerBin/(nUniquePreludeFrames*stimsetTiming.dotFramesPerCycle);
        nStepFrames = stimsetTiming.framesPerStep/nCyclesPerStepPrelude;
        if (nUniquePreludeFrames > 0)
            UniqueFramesPerStepPrelude = nUniquePreludeFrames/stimsetTiming.nPreludeBins;
        
            for s = 1:stimsetTiming.nPreludeBins
                stepFrames = zeros(nStepFrames, 1);
                stepFrames(1:stimsetTiming.dotFramesPerCycle:end) = nUniqueFrames + ((s - 1)*UniqueFramesPerStepPrelude + 1:s*UniqueFramesPerStepPrelude);
                rImSeqPrelude = cat(1, rImSeqPrelude, repmat(stepFrames, [nCyclesPerStepPrelude 1]));
            end
            
        %% for smooth prelude to stumulus to postlude transition:
        end
        rImSeq = uint32(cat(1, rImSeqPrelude, rImSeq, rImSeqPrelude));
        rIms = cat(4, rIms, rPrelude);
        rIms = uint8(255*rIms);
        
        isSuccess = true;
        output = { isSuccess, rIms, cast( rImSeq, 'int32') }; % "Images" (single) and "Image Sequence" (Int32)
        %clear rIms
        assignin( 'base', 'output', output ) 
    end
    %%%%%%%%% STIMULUS GENERATION PART
    
end
