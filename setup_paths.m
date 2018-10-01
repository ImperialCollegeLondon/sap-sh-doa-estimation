% setup_paths.m adds the required paths and does some checking

% we have an external dependecy on voicebox, which is available from
% http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
PATH_TO_VOICEBOX = fullfile('..','voicebox');

required_subdirectories = {'src','example_audio'};
for ii = 1:length(required_subdirectories)
    if ~isequal(exist(required_subdirectories{ii},'dir'),7)
        error('Missing %s subdirectory',required_subdirectories{ii});
    end
end

required_paths = {'src',PATH_TO_VOICEBOX};
for ii = 1:length(required_paths)
    if ~isequal(exist(required_paths{ii},'dir'),7)
        error('Missing %s path',required_subdirectories{ii});
    end
    addpath(required_paths{ii})
end





