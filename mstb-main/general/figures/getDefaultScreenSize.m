function [ fullscreen, doResize, def] = getDefaultScreenSize()
%getDefaultScreenSize - uses the defaults mat file in the root directory to
%create/store the screensize. This needs to be unique for everyuser, so
%must also save that information too...

% Everything is now done in normalized units.
doResize = 1;
set(0,'Units','normalized');
fullscreen = get(0,'ScreenSize'); 
def        = open('msidefaults.mat');
return;
% Determine who is the user of the toolbox, not for anything Orwellian, but
% just to ensure that default parameters are distinct for each computer.
dnam = deblank(evalc('!whoami'));
user = charPurge(dnam,'-');

% Open defaults where it might already exist
def = open('defaults.mat');

% Find if the user in defaults matches the current user...
doResize = 1;
if isfield(def.defaults,'user')    
    % See if they match?
    usermatch = strcmp(def.defaults.user,user);
    if usermatch == 1
        % If they match, then get the screensize, either from default or by
        % having to generate it the first time
        if isfield(def.defaults,'screensize')    
            fullscreen = def.defaults.screensize;
            doResize = 0;
        else
            fullscreen = get(0,'ScreenSize');
        end
    else
        % Users don't match, so get screensize now
        fullscreen = get(0,'ScreenSize');
    end
else
    % So no user in the defaults
    def.defaults.user = user;
    
    % Now get the size
    fullscreen = get(0,'ScreenSize');    
end

end