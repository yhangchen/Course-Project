rand("seed", 1);
name = "fragment2"; % file name
format = "flac";     % file format

[audio_mat, Fs] = audioread("audio\" + name + '.' + format);
% read file

for i = 1: 2
    writetable(table(audio_mat(:, i)), "data\" + name + num2str(i) + ".csv");
end
% write csv

audio_mat1 = audio_mat + randn(size(audio_mat)) * 0.01;
% add noise

for i = 1: 2
    writetable(table(audio_mat1(:, i)),...
        "data\" + name + "_noise" + num2str(i) + ".csv");
end
% write csv