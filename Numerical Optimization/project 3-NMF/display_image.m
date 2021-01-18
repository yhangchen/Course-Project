function display_image(size,fea)
% from http://www.cad.zju.edu.cn/home/dengcai/Data/Yale/images.html
%===========================================
faceW = size;
faceH = size;
numPerLine = 8;
ShowLine = 2;

Y = zeros(faceH*ShowLine,faceW*numPerLine);
for i=0:ShowLine-1
   for j=0:numPerLine-1
     Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(fea(i*numPerLine+j+1,:),[faceH,faceW]);
   end
end
imagesc(Y);colormap(gray);
%===========================================