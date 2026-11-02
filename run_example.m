% Run the estimation example with the supplied data.
cd(fileparts(mfilename('fullpath')));
load('DATA_r.mat','DATA_r');
DATA = DATA_r;
clear DATA_r;
prevalence_example;
