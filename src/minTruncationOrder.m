function[sphHarmOrderRequired] = minTruncationOrder(mic,fs,c,error_threshold_db)
% use freqMax and array geometry to determine the required
% spherical harmonic order
% k = 2*pi*f/c
% calculate truncation error at highest frequency and with sensor
% furthest from the origin based on Jin2014

if nargin < 4
    error_threshold_db = -80;
end

if nargin < 3 || isempty(c)
    c = soundspeed;
end

% Fixed upper limit on spherical order to consider
% set ridiculuously high, so error is relative to this
SPHERICAL_ORDER_UPPER_LIMIT = 500;

% Nyquist frequency sets upper frequency where worst case limit is
% determined
freqMax = fs/2;



b = modeStrength('rigid' ....
    ,mic.a ...
    ,mic.a ...
    ,freqMax,SPHERICAL_ORDER_UPPER_LIMIT,c);

% value of SPHERICAL_ORDER_UPPER_LIMIT should ensure that mode
% stregnth gets so small the values are nans
i_first_invalid = find(isnan(b),1,'first'); % values get so small they become nans
if isempty(i_first_invalid)
    % sanity check that we have gone high enough
    error('None of the mode strenghts were nan - try increasing the value of SPHERICAL_ORDER_UPPER_LIMIT in the code')
end
% remove the nan values
b(i_first_invalid:end) = [];


tot = cumsum(b .* (2*(0:length(b)-1)+1));
truncationError = abs(1-tot./tot(end));

% find order which meets predefined threshold ('1+..' because 0 order is index 1)
sphHarmOrderRequired = -1+find(truncationError<10^(error_threshold_db/20),1,'first');

if nargout==0
    %No outputs so plot the data
    figure
    plot(0:length(b)-1,20*log10(truncationError));
    hold all
    yvals = get(gca,'ylim');
    plot(sphHarmOrderRequired([1 1]),yvals,'k:')
    plot([0 length(b)-1],...
        20*log10(truncationError(1+[sphHarmOrderRequired sphHarmOrderRequired])),...
        ':k')
    plot([0 length(b)-1],...
        error_threshold_db([1 1]),...
        'r','linewidth', 2)
    xlabel('Truncation order')
    ylabel('Truncation error [dB]')
    title('Error due to truncating spherical harmonic series')
end