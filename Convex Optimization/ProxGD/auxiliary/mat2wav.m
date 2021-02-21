rand("seed", 1);
name = "fragment1"; % file name
format = "flac";     % file format
Fs = 44100;

for i = 1: 2
    audio_mat(:, i) = importdata("data\" + name + num2str(i) + ".csv");
end
% read csv

audiowrite("audio\" + name + "_rewrite" + ".wav", audio_mat, Fs);