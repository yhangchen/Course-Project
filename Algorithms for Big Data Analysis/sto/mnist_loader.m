function data = mnist_loader
data = [];
if ~exist('data\mnist', 'dir')
  mkdir('data\mnist') ;
end
names = {'train-images-idx3-ubyte', 'train-labels-idx1-ubyte', ...
         't10k-images-idx3-ubyte', 't10k-labels-idx1-ubyte'} ;
for i=1:4
  if ~exist(fullfile('data\mnist', names{i}), 'file')
    url = sprintf('http://yann.lecun.com/exdb/mnist/%s.gz',names{i}) ;
    fprintf('downloading %s\n', url) ;
    gunzip(url, 'data\mnist') ;
  end
end
data.train_img = loadMNISTImages('data\mnist\train-images-idx3-ubyte');
data.train_lbl = loadMNISTLabels('data\mnist\train-labels-idx1-ubyte');
data.test_img = loadMNISTImages('data\mnist\t10k-images-idx3-ubyte');
data.test_lbl = loadMNISTLabels('data\mnist\t10k-labels-idx1-ubyte');
end

% The following two functions are from 
% https://github.com/amaas/stanford_dl_ex/tree/master/common

function images = loadMNISTImages(filename)
fp = fopen(filename, 'rb');
assert(fp ~= -1, ['Could not open ', filename, '']);
magic = fread(fp, 1, 'int32', 0, 'ieee-be');
assert(magic == 2051, ['Bad magic number in ', filename, '']);
numImages = fread(fp, 1, 'int32', 0, 'ieee-be');
numRows = fread(fp, 1, 'int32', 0, 'ieee-be');
numCols = fread(fp, 1, 'int32', 0, 'ieee-be');
images = fread(fp, inf, 'unsigned char');
images = reshape(images, numCols, numRows, numImages);
images = permute(images,[2 1 3]);
fclose(fp);
end

function labels = loadMNISTLabels(filename)
fp = fopen(filename, 'rb');
assert(fp ~= -1, ['Could not open ', filename, '']);
magic = fread(fp, 1, 'int32', 0, 'ieee-be');
assert(magic == 2049, ['Bad magic number in ', filename, '']);
numLabels = fread(fp, 1, 'int32', 0, 'ieee-be');
labels = fread(fp, inf, 'unsigned char');
assert(size(labels,1) == numLabels, 'Mismatch in label count');
fclose(fp);
end

