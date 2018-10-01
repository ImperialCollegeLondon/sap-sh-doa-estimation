function[X,pm,t,f,out_tail] = stft(x,w,inc,n_fft,fs,in_tail)
%STFT takes 2D multichannel signal matrix and returns 3D short time fourier
%transform matrix at n_freq=fix(1+n_fft/2) frequencies
%
%[X,pm,t,f] = stft(x,w,inc,n_fft,fs)
%
% Inputs
%       x  input signal [n_samples n_channels]
%      fs  sample rate of input signal
%       w  frame length or window of frame length
%     inc  frame increment in samples
%
% Outputs
%       X  STFT matrix [n_freq, n_channels, n_frames]
%       t  centre time of each frame
%       f  frequency frequency of each bin
%
%

% %% dummy data for testing
% x = bsxfun(@times,[repmat([1:20]',1,4)],[1 10 100 1000])
% fs = 1;
% w = ones(8,1);
% inc = 4;

FORCE_POW_2 = 1;
% *** TODO generalise to have arbitrary window length
if ~isempty(n_fft) && ~isscalar(n_fft)
    error('n_fft must be a scalar')
end


if nargin < 4 || isempty(n_fft)
    n_fft = length(w);
    if n_fft == 1
        n_fft = w;
    end
end



if nargin < 6
    block_based = 0;
else
    block_based = 1;
end
    
if isscalar(w)
    if w/inc ~= 2
        error('Overlap factor is not two so full window must be specified')
    end
    w = sqrt(hamming(w,'periodic'));
    w=w/sqrt(sum(w(1:inc:length(w)).^2)); %normalise window for overlap factor of 2
end

n_w = length(w);

if FORCE_POW_2 && 2^nextpow2(n_fft)~=n_fft
    %error('Window must be power of 2')
    warning('Making n_fft a power of 2')
    n_fft = 2^nextpow2(n_fft);
end

if n_fft > n_w % will apply padding to beginning and end    
    fft_pad = n_fft-n_w;
    fft_pre_pad = round(fft_pad/2);
    fft_post_pad = fft_pad-fft_pre_pad;
    else
    fft_pre_pad = 0;
    fft_post_pad = 0;
end

[n_samples,n_chans] = size(x);

if ~block_based
    %no tail to conisder so pad with zeros
    pre_pad_len = n_w-inc;                                % last inc samples of first frame will be non-zero
    post_pad_len = inc-mod(n_samples,inc) + pre_pad_len;    % fill final block and then pad
    x = [zeros(pre_pad_len,n_chans);...
        x;...
        zeros(post_pad_len,n_chans)];
else
    %prepend with saved tail
    pre_pad_len=0;
    post_pad_len=0;
    x = [in_tail; x];

end

len_x = size(x,1);

% enframe signal
fr_st = 1:inc:len_x-(n_w-1);
n_frames = length(fr_st);
fr_idc = bsxfun(@plus,fr_st,[0:n_w-1]');

if block_based
    out_tail = x(fr_st(end)+inc:end,:); %tail starts with first incomplete frame
end

y = permute(reshape(x(fr_idc,:),[n_w,n_frames,n_chans]),[1,3,2]);
%pm.y = y;

% apply window
y = y .* repmat(w(:),[1, n_chans, n_frames]);

if fft_pre_pad
    y = [zeros(fft_pre_pad,n_chans,n_frames); y; zeros(fft_post_pad,n_chans,n_frames)];
end

% do fft
X = rfft(y,n_fft,1);

% bundle up parameters required for istft
pm.inc = inc;
pm.n_fft = n_fft;
pm.w = w;
pm.fr_st = fr_st;
pm.fr_idc = fr_idc;
pm.pre_pad_len = pre_pad_len;
pm.post_pad_len = post_pad_len;
pm.len_x = len_x;
pm.n_samples = n_samples;
pm.fft_pre_pad = fft_pre_pad;
pm.fft_post_pad = fft_post_pad;


% calculate time and frequency scales
if nargin < 4 || isempty(fs)
    fs = 1;
end
t = (fr_st-1 - (pre_pad_len) + n_fft/2) ./fs;
f = [0:n_fft/2] * fs/n_fft;
pm.t = t;
pm.f = f;
pm.fs = fs;

