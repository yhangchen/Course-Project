function out = test4(tol)
pic = imread('fig_pic.jpg');
imshow(pic(:,:,:))
pic = im2double(pic);
[nrow, ncol] = size(pic(:,:,1));
pic = [pic(:,:,1);pic(:,:,2);pic(:,:,3)];
rank(pic)
sig=0.01;
prop_obs=.2;
R=sprand(nrow,ncol,prop_obs);
[ii, jj]= find(R);
R = sparse(ii, jj, 1, nrow,ncol);
R = [R;R;R];
xobs=(pic+sig*randn(nrow*3,ncol)).*R;
opics = zeros(nrow,ncol,3);
opics(:,:,1) = xobs(1:nrow,:);
opics(:,:,2) = xobs(nrow+1:2*nrow,:);
opics(:,:,3) = xobs(2*nrow+1:3*nrow,:);
imshow(opics)
fprintf('Starting SVT\n')
out = svt(R,xobs,0,tol,[0,255]);
fprintf('SVT end\n\n')

X = out.X;
test_err = norm(X-pic,'fro')/norm(pic,'fro')
pics = zeros(nrow,ncol,3);
pics(:,:,1) = X(1:nrow,:);
pics(:,:,2) = X(nrow+1:2*nrow,:);
pics(:,:,3) = X(2*nrow+1:3*nrow,:);
imshow(pics)

