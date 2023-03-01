%----------------------------Face_LineJudgement---------------------------%
%AUTHOR INFORMATION: Jake J. Son, Nate M. Petro
%Institute for Human Neuroscience, Boys Town National Research Hospital
%INPUTS: Participant Number

%-------------------------------------------------------------------------%
%                             Trigger Legend                              %
%-------------------------------------------------------------------------%
% fixation = 60

% small angle difference, both top = 41
% small angle difference, both bottom = 42

% large angle difference, both top = 81
% large angle difference, both bottom = 82

%-------------------------------------------------------------------------%
%                             General Setup                               %
%-------------------------------------------------------------------------%

sca;
clear all;
clear vars;
% com1 = serialport('COM1',19200);

rand('state', sum(100*clock));  %set starting state of clock randomly each time - keeps randomization calls from repeating%
Screen('Preference', 'Verbosity', 3);           %-  was 3
Screen('Preference', 'SkipSyncTests', 2);       % Set the verbosity and extent of initialization tests and debugging%
Screen('Preference', 'VisualDebugLevel', 1);    

%-------------------------------------------------------------------------%
%                             Color Scheme                                %
%-------------------------------------------------------------------------%
gray = [127 127 127]; white = [255 255 255]; black = [0 0 0];
green = [0 255 0]; red = [255 0 0];

bgcolor = black; textcolor = gray;

alpha=1;

KbName('UnifyKeyNames');
escKey = KbName('ESCAPE');

inputs = {'Participant number'};
defaults = {'0'};	%prompt for experiment parameters, read inputs into variables%
answer = inputdlg(inputs, 'Face_LineJudgement', 2, defaults);
[PartID] = deal(answer{:});

directory = fileparts(mfilename('fullpath'));

% DefaultName = strcat('Speedball_', PartID);
% PathName=strcat(directory,'\output');
% FileName=strcat(DefaultName,'.xls');
% outputfile = fopen(fullfile(PathName,FileName),'w+'); 						  
% fprintf(outputfile, 'PartId\t Trial Number\t Trigger\t Jitter (s)\n');

[win, winDimen] = Screen('OpenWindow', max(Screen('Screens')));
Screen('FillRect', win, bgcolor);
center = [winDimen(3)/2, winDimen(4)/2];
Screen('Flip',win);

% ------------------------------------------------------------------------%
%                        Initialize task parameters                       %
% ------------------------------------------------------------------------%

FixDur = 2.0;
TrialNum = 200;

% ------------------------------------------------------------------------%
%                        Initialize line positions                        %
% ------------------------------------------------------------------------%

%Top = Left; Bottom = Right;
straightLineTop_point1x    = (1/8)*winDimen(3)*-1 ;
straightLineTop_point1y    = (1/12)*winDimen(4)*-1 ;
straightLineTop_point2x    = (1/8)*winDimen(3)*-1  ;
straightLineTop_point2y    = (1/12)*winDimen(4) ;

straightLineBottom_point1x = (1/8)*winDimen(3) ;
straightLineBottom_point1y = (1/12)*winDimen(4)*-1 ;
straightLineBottom_point2x = (1/8)*winDimen(3)  ;
straightLineBottom_point2y = (1/12)*winDimen(4) ;

HideCursor();

%-------------------------------------------------------------------------%
%                            ProPixx Dot Parameters                       %
%-------------------------------------------------------------------------%

dotColorBlack = black;
dotColorWhite = white;
dotSizePix = 5; 
    
%-------------------------------------------------------------------------%
%                              Load Task List                             %
%-------------------------------------------------------------------------%  

table = fullfile(directory, 'stimList_FaceLine_Practice.txt'); 
opts = detectImportOptions(table); %setup table import from selected
trial_index = readtable(table,opts);
clear opts;

Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

Priority(1);

%-------------------------------------------------------------------------%
%                          Fixation Parameters                            %
%-------------------------------------------------------------------------%

fixCrossDimPix = 0.02*winDimen(4); %sets the fixation cross arm size (this is a ratio of the window dimensions--3% of the y window dimension)
fixWidth = round(.0033*winDimen(4));
fixLineCoords = [-fixCrossDimPix fixCrossDimPix 0 0; 0 0 -fixCrossDimPix fixCrossDimPix];   %set fixation stimulus parameters relative to zero based on pixel/ratio size (for the x plane, how many pixels left/right the cross will take up; for the y plane, how many pixels up/down the cross will take up)
JitterLengths = (randi([(FixDur*1000)-300,(FixDur*1000)+300],1,TrialNum))./1000; %sets a random jitter length +/- 200 ms from our fixation duration length (3.75 ms)--should be 10% of your trial length

%-------------------------------------------------------------------------%
%                             Trial Loop                                  %                            
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%                             Fixation                                    %
%-------------------------------------------------------------------------%

Screen('DrawLines',win,fixLineCoords,3,textcolor,center);
Screen('DrawDots', win, [0 0], dotSizePix, dotColorBlack, [], 2);
Screen('DrawingFinished', win);  					
Screen('Flip',win);
% write(com1,60,"uint16");
WaitSecs(FixDur);

for i=1:TrialNum
    
        [KeyIsDown,~,KeyID] = KbCheck;
        if KeyIsDown
            if KeyID(escKey)
                sca;
                clearvars;
                % clear com1;
                error('Experiment aborted by user!')
               
            end
        end

    if isnumeric(trial_index.Trigger(i))==1

%-------------------------------------------------------------------------%
%                           Set degree to rotate BOTH LINES               %
%-------------------------------------------------------------------------%               

        if strcmp(trial_index.Parallel(i),"P")==1 & strcmp(trial_index.Line(i),"C")
            [xf_top, yf_top] = lrotate(...
            [straightLineTop_point1x straightLineTop_point2x], ... % x coordinates of the 2 points in the line
            [straightLineTop_point1y straightLineTop_point2y], ... % y coordinates of the 2 points in the line
            trial_index.Origin(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineTop_point1x + [abs(straightLineTop_point1x - straightLineTop_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineTop_point1y + [abs(straightLineTop_point1y - straightLineTop_point2y) / 2] ...  % y coordinate about which to rotate the line
            );

            [xf_bot, yf_bot] = lrotate(...
            [straightLineBottom_point1x straightLineBottom_point2x], ... % x coordinates of the 2 points in the line
            [straightLineBottom_point1y straightLineBottom_point2y], ... % y coordinates of the 2 points in the line
            trial_index.Origin(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineBottom_point1x + [abs(straightLineBottom_point1x - straightLineBottom_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineBottom_point1y + [abs(straightLineBottom_point1y - straightLineBottom_point2y) / 2] ...  % y coordinate about which to rotate the line
            );

        elseif strcmp(trial_index.Parallel(i),"P")==1 & strcmp(trial_index.Line(i),"R")
            [xf_top, yf_top] = lrotate(...
            [straightLineTop_point1x straightLineTop_point2x], ... % x coordinates of the 2 points in the line
            [straightLineTop_point1y straightLineTop_point2y], ... % y coordinates of the 2 points in the line
            -trial_index.Origin(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineTop_point1x + [abs(straightLineTop_point1x - straightLineTop_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineTop_point1y + [abs(straightLineTop_point1y - straightLineTop_point2y) / 2] ...  % y coordinate about which to rotate the line
            );

            [xf_bot, yf_bot] = lrotate(...
            [straightLineBottom_point1x straightLineBottom_point2x], ... % x coordinates of the 2 points in the line
            [straightLineBottom_point1y straightLineBottom_point2y], ... % y coordinates of the 2 points in the line
            -trial_index.Origin(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineBottom_point1x + [abs(straightLineBottom_point1x - straightLineBottom_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineBottom_point1y + [abs(straightLineBottom_point1y - straightLineBottom_point2y) / 2] ...  % y coordinate about which to rotate the line
            );


        elseif strcmp(trial_index.Parallel(i),"N") & strcmp(trial_index.TargetSide(i),"L") == 1 & strcmp(trial_index.Line(i),"C") %choose whether to rotate the top or bottom line????
            [xf_top, yf_top] = lrotate(...
            [straightLineTop_point1x straightLineTop_point2x], ... % x coordinates of the 2 points in the line
            [straightLineTop_point1y straightLineTop_point2y], ... % y coordinates of the 2 points in the line
            trial_index.Origin(i) + trial_index.Tilt(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineTop_point1x + [abs(straightLineTop_point1x - straightLineTop_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineTop_point1y + [abs(straightLineTop_point1y - straightLineTop_point2y) / 2] ...  % y coordinate about which to rotate the line
            );

            [xf_bot, yf_bot] = lrotate(...
            [straightLineBottom_point1x straightLineBottom_point2x], ... % x coordinates of the 2 points in the line
            [straightLineBottom_point1y straightLineBottom_point2y], ... % y coordinates of the 2 points in the line
            trial_index.Origin(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineBottom_point1x + [abs(straightLineBottom_point1x - straightLineBottom_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineBottom_point1y + [abs(straightLineBottom_point1y - straightLineBottom_point2y) / 2] ...  % y coordinate about which to rotate the line
            );
        elseif strcmp(trial_index.Parallel(i),"N") & strcmp(trial_index.TargetSide(i),"L") == 1 & strcmp(trial_index.Line(i),"R") %choose whether to rotate the top or bottom line????
            [xf_top, yf_top] = lrotate(...
            [straightLineTop_point1x straightLineTop_point2x], ... % x coordinates of the 2 points in the line
            [straightLineTop_point1y straightLineTop_point2y], ... % y coordinates of the 2 points in the line
            -trial_index.Origin(i) - trial_index.Tilt(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineTop_point1x + [abs(straightLineTop_point1x - straightLineTop_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineTop_point1y + [abs(straightLineTop_point1y - straightLineTop_point2y) / 2] ...  % y coordinate about which to rotate the line
            );

            [xf_bot, yf_bot] = lrotate(...
            [straightLineBottom_point1x straightLineBottom_point2x], ... % x coordinates of the 2 points in the line
            [straightLineBottom_point1y straightLineBottom_point2y], ... % y coordinates of the 2 points in the line
            -trial_index.Origin(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineBottom_point1x + [abs(straightLineBottom_point1x - straightLineBottom_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineBottom_point1y + [abs(straightLineBottom_point1y - straightLineBottom_point2y) / 2] ...  % y coordinate about which to rotate the line
            );
        elseif strcmp(trial_index.Parallel(i),"N") & strcmp(trial_index.TargetSide(i),"R") == 1 & strcmp(trial_index.Line(i),"C") %choose whether to rotate the top or bottom line????
            [xf_top, yf_top] = lrotate(...
            [straightLineTop_point1x straightLineTop_point2x], ... % x coordinates of the 2 points in the line
            [straightLineTop_point1y straightLineTop_point2y], ... % y coordinates of the 2 points in the line
            trial_index.Origin(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineTop_point1x + [abs(straightLineTop_point1x - straightLineTop_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineTop_point1y + [abs(straightLineTop_point1y - straightLineTop_point2y) / 2] ...  % y coordinate about which to rotate the line
            );

            [xf_bot, yf_bot] = lrotate(...
            [straightLineBottom_point1x straightLineBottom_point2x], ... % x coordinates of the 2 points in the line
            [straightLineBottom_point1y straightLineBottom_point2y], ... % y coordinates of the 2 points in the line
            trial_index.Origin(i) + trial_index.Tilt(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineBottom_point1x + [abs(straightLineBottom_point1x - straightLineBottom_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineBottom_point1y + [abs(straightLineBottom_point1y - straightLineBottom_point2y) / 2] ...  % y coordinate about which to rotate the line
            );

        elseif strcmp(trial_index.Parallel(i),"N") & strcmp(trial_index.TargetSide(i),"R") == 1 & strcmp(trial_index.Line(i),"R") %choose whether to rotate the top or bottom line????
            [xf_top, yf_top] = lrotate(...
            [straightLineTop_point1x straightLineTop_point2x], ... % x coordinates of the 2 points in the line
            [straightLineTop_point1y straightLineTop_point2y], ... % y coordinates of the 2 points in the line
            -trial_index.Origin(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineTop_point1x + [abs(straightLineTop_point1x - straightLineTop_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineTop_point1y + [abs(straightLineTop_point1y - straightLineTop_point2y) / 2] ...  % y coordinate about which to rotate the line
            );

            [xf_bot, yf_bot] = lrotate(...
            [straightLineBottom_point1x straightLineBottom_point2x], ... % x coordinates of the 2 points in the line
            [straightLineBottom_point1y straightLineBottom_point2y], ... % y coordinates of the 2 points in the line
            -trial_index.Origin(i) - trial_index.Tilt(i), ...                                            %angle (in degrees) to rotate the straight line
            straightLineBottom_point1x + [abs(straightLineBottom_point1x - straightLineBottom_point2x) / 2],...  % x coordinate about which to rotate the line
            straightLineBottom_point1y + [abs(straightLineBottom_point1y - straightLineBottom_point2y) / 2] ...  % y coordinate about which to rotate the line
            );

        end

        % ------------------------------
        % Initialize movement parameters
        % ------------------------------

    fps=Screen('FrameRate', win);      % frames per second
    ifi=Screen('GetFlipInterval', win);

        if fps==0
           fps=1/ifi;
        end

    LinesPointsVec(1,1:4) =  [xf_top(1),  xf_top(2), xf_bot(1),  xf_bot(2)]; %set up lines from previous loop
    LinesPointsVec(2,1:4) =  [yf_top(1),  yf_top(2), yf_bot(1),  yf_bot(2)];

    %% Experimental 

    nframes1=90; %determines duration of movement (framerate assumed to be 60Hz)
    nframes2=60;
    % Reduce nframes to increase speed of moving ball

    dx1 = (xf_top(2) - xf_top(1))/nframes1;
    dy1 = (yf_top(2) - yf_top(1))/nframes1;
    dx2 = (xf_bot(2) - xf_bot(1))/nframes2;
    dy2 = (yf_bot(2) - yf_bot(1))/nframes2;

    nframes_max = max(nframes1, nframes2);

    dxdy = [dx1 dx2; dy1 dy2];

    start_pos = 'bot'; % for testing purposes. will be included in new stimulus file

        if strcmp(start_pos, 'top')

            xymatrix = LinesPointsVec(1:2, 1:2:end); % determine which points are starting points based on up / down movement

        else

            xymatrix = LinesPointsVec(1:2, 2:2:end); % determine which points are starting points based on up / down movement
            dxdy = -dxdy;

        end

    vbl = Screen('Flip', win);

        for i = 1:nframes_max

            if (i>1 & xymatrix(1, 1) < xf_top(2) & xymatrix(1, 1) > xf_top(1) & xymatrix(1,2) < xf_bot(2) & xymatrix(1,2) > xf_bot(1))

                Screen('DrawDots', win, xymatrix, dotSizePix, dotColorWhite, center, 1);  % change 1 to 0 or 4 to draw square dots

                Screen('DrawingFinished', win); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            end

            xymatrix = xymatrix + dxdy; % move dots

            waitframes=1;

            vbl=Screen('Flip', win, vbl + (waitframes-0.5)*ifi);

        end

    %   % write(com1,trial_index.Trigger(i),"uint16");

    %%
    %------------------------------------------------------------------------%
    %                               Fixation                                 %  
    %------------------------------------------------------------------------%
        Screen('DrawLines',win,fixLineCoords,3,textcolor,center);
        Screen('DrawDots', win, [0 0], dotSizePix, dotColorBlack, [], 2);
        Screen('DrawingFinished', win);  					
        Screen('Flip',win);
        % write(com1,60,"uint16");
        WaitSecs(JitterLengths(1,i));



    end 
 
end

clearvars;
Priority(0);
sca;      