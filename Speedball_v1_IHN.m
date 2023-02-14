%----------------------------Face_LineJudgement---------------------------%
%AUTHOR INFORMATION: Jake J. Son, Nate M. Petro
%Institute for Human Neuroscience, Boys Town National Research Hospital
%INPUTS: Participant Number

%-------------------------------------------------------------------------%
%                                  Trigger Legend                          %
%-------------------------------------------------------------------------%
% fixation = 60

% small angle difference, both top = 41
% small angle difference, both bottom = 42
% small angle difference, left top right bottom = 43
% small angle difference, left bottom right top = 44

% large angle difference, both top = 81
% large angle difference, both bottom = 82
% large angle difference, left top right bottom = 83
% large angle difference, left bottom right top = 84

%------------------------------------------------------------------------%
%                             General Setup                              %
%------------------------------------------------------------------------%
sca;
clear all;
clear vars;
% com1 = serialport('COM1',19200);

rand('state', sum(100*clock));  %set starting state of clock randomly each time - keeps randomization calls from repeating%
Screen('Preference', 'Verbosity', 3);           %-  was 3
Screen('Preference', 'SkipSyncTests', 2);       % Set the verbosity and extent of initialization tests and debugging%
Screen('Preference', 'VisualDebugLevel', 1);    

%------------------------------------------------------------------------%
%                             Color Scheme                               %
%------------------------------------------------------------------------%
gray = [127 127 127]; white = [255 255 255]; black = [0 0 0]; green = [0 255 0]; red = [255 0 0]; %set color indices and bgcolor%
bgcolor = black; textcolor = gray;
alpha=1;

KbName('UnifyKeyNames');
escKey = KbName('ESCAPE');
%     
%     inputs = {'Participant number'};
%     defaults = {'0'};	%prompt for experiment parameters, read inputs into variables%
%     answer = inputdlg(inputs, 'Face_LineJudgement', 2, defaults);
%     [PartID] = deal(answer{:});

directory = fileparts(mfilename('fullpath'));

%     DefaultName = strcat('Face_LineJudgement_', PartID);
%     PathName=strcat(directory,'\output');
%     FileName=strcat(DefaultName,'.xls');
%     outputfile = fopen(fullfile(PathName,FileName),'w+'); 						  
%     fprintf(outputfile, 'PartId\t Trial Number\t Trigger\t Jitter (s)\n');
%     
[win, winDimen] = Screen('OpenWindow', max(Screen('Screens')));
Screen('FillRect', win, bgcolor);
center = [winDimen(3)/2, winDimen(4)/2];	%set screen parameters and important stimulus positions%
Screen('Flip',win);

FixDur = 2.0;
ProbeDur = 1.5;
TrialNum = 200;
resizePercent = .5;

%Top = Left; Bottom = Right;
straightLineTop_point1x    = (1/8)*winDimen(3)*-1 ;
straightLineTop_point1y    = (1/12)*winDimen(4)*-1 ;
straightLineTop_point2x    = (1/8)*winDimen(3)*-1  ;
straightLineTop_point2y    = (1/12)*winDimen(4) ;

straightLineBottom_point1x = (1/8)*winDimen(3) ;
straightLineBottom_point1y = (1/12)*winDimen(4)*-1 ;
straightLineBottom_point2x = (1/8)*winDimen(3)  ;
straightLineBottom_point2y = (1/12)*winDimen(4) ;

colvect = ones(3,4) .* gray(1); % make both lines gray (will be called later)

HideCursor();

%-------------------------------------------------------------------------%
%                            ProPixx Dot Parameters                       %
%-------------------------------------------------------------------------%
    dotColorBlack = black;
    dotColorWhite = white;
    dotSizePix = 2; 
    
%-------------------------------------------------------------------------%
%                                Load Task List                           %
%-------------------------------------------------------------------------%  
 
    table = fullfile(directory, 'stimList_FaceLine_Practice.txt');  
    opts = detectImportOptions(table); %setup table import from selected
    trial_index = readtable(table,opts);
    clear opts;
    list = dir(fullfile(directory,'\stims','*.tif')); %---assigns the stimulus directory
    filenames = {list.name}; %---makes the list of image file names a string
    
    Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    Priority(1);   
    
    
%     trial_index.

%-------------------------------------------------------------------------%
%                                     Fixation Parameters                                           %
%-------------------------------------------------------------------------%
    fixCrossDimPix = 0.02*winDimen(4); %sets the fixation cross arm size (this is a ratio of the window dimensions--3% of the y window dimension)
    fixWidth = round(.0033*winDimen(4));
    fixLineCoords = [-fixCrossDimPix fixCrossDimPix 0 0; 0 0 -fixCrossDimPix fixCrossDimPix];   %set fixation stimulus parameters relative to zero based on pixel/ratio size (for the x plane, how many pixels left/right the cross will take up; for the y plane, how many pixels up/down the cross will take up)
    JitterLengths = (randi([(FixDur*1000)-300,(FixDur*1000)+300],1,TrialNum))./1000; %sets a random jitter length +/- 200 ms from our fixation duration length (3.75 ms)--should be 10% of your trial length

    face_stims = dir(fullfile(directory,'\stims','*.tif'));

%-------------------------------------------------------------------------%
                               %Trial Loop%                            
%-------------------------------------------------------------------------%

%------------------------------------------------------------------------%
%                             Fixation                                   %
%------------------------------------------------------------------------%
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
%         if i == round(.5*TrialNum)%%% 15s break after 1/2 of trials complete %%%
%            Screen('FillRect',win,gray);
%            Screen('TextSize',win,24);
%            DrawFormattedText(win, 'Great Job!\n Please remain still and feel free to rest your eyes.\n The task will resume in 15 seconds.','center','center',black); 
%            Screen('Flip',win);
%            WaitSecs(15);
%            Screen('FillRect',win,black);
% 
%            
%         end
        
        
if isnumeric(trial_index.Trigger(i))==1 
                
%-------------------------------------------------------------------------%
%                           Set degree to rotate BOTH LINES
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



LinesPointsVec(1,1:4) =  [xf_top(1),  xf_top(2), xf_bot(1),  xf_bot(2)]; %set up lines from previous loop
LinesPointsVec(2,1:4) =  [yf_top(1),  yf_top(2), yf_bot(1),  yf_bot(2)];


  Face = trial_index.FaceFile(i);   %find face file for this trial
  base_path = list(1).folder;
  image_path_face = fullfile(base_path,Face);
  trial_image_face = imread(char(image_path_face));
  trial_image_face = imresize(trial_image_face,resizePercent);
  trial_image_texture_face = Screen('MakeTexture',win, trial_image_face);
  FaceStim_rect=Screen('Rect', trial_image_texture_face);
  Screen('DrawTexture',win,trial_image_texture_face, [],CenterRect(FaceStim_rect, winDimen)); %draw face texture on screen
  Screen('DrawLines', win, LinesPointsVec , 5, colvect, center, 1); %draw lines from previous loop setup
  Screen('DrawDots', win, [0 0], dotSizePix, dotColorWhite, [], 2); %propixx dot
  Screen('DrawingFinished',win);
  Screen('Flip',win);
  % write(com1,trial_index.Trigger(i),"uint16");
  WaitSecs(ProbeDur);


%------------------------------------------------------------------------%
%                               Fixation         %  
%------------------------------------------------------------------------%
    Screen('DrawLines',win,fixLineCoords,3,textcolor,center);
    Screen('DrawDots', win, [0 0], dotSizePix, dotColorBlack, [], 2);
    Screen('DrawingFinished', win);  					
    Screen('Flip',win);
    % write(com1,60,"uint16");
    WaitSecs(JitterLengths(1,i));



end %trial loop
  
%     fprintf(outputfile, '%s\t %d\t %d\t %.5f\n', PartID, i, trial_index.Trigger(i), JitterLengths(1,i));
 
end

% clear com1;
% fclose(outputfile);
clearvars;
Priority(0);
sca;      