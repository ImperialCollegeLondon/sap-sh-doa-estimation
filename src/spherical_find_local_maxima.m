function[inc,iaz] = spherical_find_local_maxima(X)

[ninc,naz] = size(X);

buf = zeros(ninc+2,naz+2);
buf(2:end-1,2:end-1) = X;
buf(2:end-1,1) = X(:,end);
buf(2:end-1,end) = X(:,1);
buf(1,2:end-1) = fliplr(X(1,:));
buf(end,2:end-1) = fliplr(X(end,:));

[inc,iaz] = find( X > buf(1:end-2,2:end-1) & ... %compare to the top
                   X > buf(3:end,2:end-1) & ...  %...bottom
                   X > buf(2:end-1,1:end-2) & ... %...left
                   X > buf(2:end-1,3:end) ); %...right