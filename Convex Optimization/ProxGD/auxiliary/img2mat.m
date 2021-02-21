rand("seed", 1);
img = imread("Lenna.png");

for i = 1: 3
    writetable(table(img(:, :, i)), "datasets\img" + num2str(i) + ".csv");
end

img1 = max(min(0,double(img) + randn(size(img)) * 255 * 0.1),255);
for i = 1: 3
    writetable(table(img1(:, :, i)), "datasets\img_noise" + num2str(i) + ".csv");
end