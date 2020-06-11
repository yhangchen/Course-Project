function data = covertype_loader
data = [];
if ~exist('data\covertype', 'dir')
  mkdir('data\covertype') ;
end
if ~exist('data\covertype\covtype.libsvm.binary', 'file')
    url = sprintf('https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/covtype.libsvm.binary.bz2') ;
    fprintf('Please download file from %s,\n and uncompress it to data\\covertype \n', url);
    % matlab cannot uncompress with .bz2 file.
    error('No data');
end
[data.label, data.data] = libsvmread('data\covertype\covtype.libsvm.binary');
end