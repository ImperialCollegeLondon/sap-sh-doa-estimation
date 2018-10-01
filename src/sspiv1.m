function[pivs, erank] = sspiv1(Ucell,Scell)

%Assume that Ucell and Scell are arranged as [nFreq, nFrames]
[nFreq, nFrames] = size(Ucell);
if ~isequal(size(Ucell), size(Scell))
    error('Ucell and Scell must have the same dimensions');
end

sigspace = NaN(nFreq,4,nFrames);
erank = NaN(nFreq,nFrames);

for i_freq = 1:nFreq
    for i_frame = 1:nFrames
        if ~isempty(Ucell{i_freq,i_frame})     
            erank(i_freq,i_frame) = Scell{i_freq,i_frame}(1)/Scell{i_freq,i_frame}(2);
            sigspace(i_freq,:,i_frame) = Ucell{i_freq,i_frame}(1:4,1);
        end
    end
end
pivs = piv(sigspace);
