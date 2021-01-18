function Results = test(r)

addpath('linesearch')
addpath('datasets\Yale')
addpath('datasets\ORL')
addpath('datasets\bbcsport')
addpath('datasets\mining')

% default r 10.
datas = {};
datas{1} = rand(2048,r)*rand(r,1024);

data_in = load('Yale_32x32.mat');
datas{2} = data_in.fea;
data_in = load('Yale_64x64.mat');
datas{3} = data_in.fea;
data_in = load('ORL_32x32.mat');
datas{4} = data_in.fea;
data_in = load('ORL_64x64.mat');
datas{5} = data_in.fea;

datas{6} = mm_to_msm('bbcsport.mtx');
data_in = load('tr11.mat');
datas{7} = data_in.V;



clc

fprintf('Dataset loaded');

Results = [];



for data_ind = 1:7
    V = datas{data_ind};
    if data_ind == 2 || data_ind == 3
        r = 15;
    elseif data_ind == 4 || data_ind == 5
        r = 25;
    elseif data_ind == 7
        r = 10;
    end
    
    opts.tol_grad = 1e-4;
%     opts.metric = 'fro';
%     [x1, infos] = nmf_mu(V, r, opts);
%     results(1,1) = infos.epoch; results(1,2) = infos.time;results(1,3) = infos.rel_cost;results(1,4) = infos.rel_projnorm;
% 
%     opts.alg = 'pgd';
%     [x2, infos] = nmf_pgd(V, r, opts);
%     results(2,1) = infos.epoch; results(2,2) = infos.time;results(2,3) = infos.rel_cost;results(2,4) = infos.rel_projnorm;
% 
% 
%     opts.alg = 'direct_pgd';
%     [x3, infos] = nmf_pgd(V, r, opts);
%     results(3,1) = infos.epoch; results(3,2) = infos.time; results(3,3) = infos.rel_cost; results(3,4) = infos.rel_projnorm;
% 
%     [x4, infos] = nmf_BB(V, r,opts);
%     results(4,1) = infos.epoch; results(4,2) = infos.time; results(4,3) = infos.rel_cost;results(4,4) = infos.rel_projnorm;
%     
%     
%     [x5, infos] = nmf_newton_inexact(V, r,opts);
%     results(5,1) = infos.epoch; results(5,2) = infos.time; results(5,3) = infos.rel_cost;results(5,4) = infos.rel_projnorm;
% 
%     opts.tol_grad = 1e-4;
%     opts.alg = 'anls_asgroup';
%     [x6, infos] = nmf_anls(V, r, opts);
%     results(6,1) = infos.epoch; results(6,2) = infos.time;  results(6,3) = infos.rel_cost;results(6,4) = infos.rel_projnorm;
% 
% 
%     opts.alg = 'anls_asgivens';
%     [x7, infos] = nmf_anls(V, r, opts);
%     results(7,1) = infos.epoch; results(7,2) = infos.time; results(7,3) = infos.rel_cost;results(7,4) = infos.rel_projnorm;
% 
%     
%     opts.alg = 'anls_bpp';
%     [x8, infos] = nmf_anls(V, r, opts);
%     results(8,1) = infos.epoch; results(8,2) = infos.time; results(8,3) = infos.rel_cost;results(8,4) = infos.rel_projnorm;
% 
    [x9, infos] = nmf_cg(V, r, opts);
    results(9,1) = infos.epoch; results(9,2) = infos.time; results(9,3) = infos.rel_cost;results(9,4) = infos.rel_projnorm;

    save(string(data_ind)+'.mat','results');
    Results = [Results;results];
end
Results