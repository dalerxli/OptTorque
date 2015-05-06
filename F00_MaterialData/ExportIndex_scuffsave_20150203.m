%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Matlab code to convert material data into scuff-em format %%%%%
%%%%%%% With an option to plot the material data 
%%%%%%% Yoonkyung Eunnie Lee 2015-03-02                         %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/LinkDropbox/Data/Matlab/elutils'); 
addpath('~/Documents/LinkDropbox/Data/Matlab/elutils'); 
SetPlotDefaults(); %%includes housekeeping / setting default styles.
    fonm = 'arial';
    lfsize = 15;
    sfsize = 12;
    %%colorlist= 'rgbkmcrgbkmc';
    extension = '.jpg';
c_light =299792458; 

basedir = '/home/eunnie12/Documents/LinkDropbox/Data/GMM/1_ExportIndex_v20150302'; 
addpath(basedir); 
%filename = 'Fe_Palik_wvnm_eps_raw.dat';
%outname = 'Fe_Palik_omega_eps_raw.dat';
dirlist = dir('wvnm_eps_raw/*_raw.dat'); 

filename='Au_Palik_wvnm_eps_raw.dat';
outname = 'Au_Palik_omega_eps_smooth_scuff2.dat';

%for (ii=1:length(dirlist))
%    filename=dirlist(ii).name; 
%    outname = filename
    if (exist(strcat('wvnm_eps_raw/',filename))==2)
        fid = fopen(strcat('wvnm_eps_raw/',filename));
        formatSpec = '%f %f %f';
        %data = fscanf(fid,format,[98 9]);
        data = textscan(fid,formatSpec,'HeaderLines',0,'CollectOutput',1,'EmptyValue',0);
        data = data{1,1};
        fclose(fid); 
    end   %%if     
%end %%for

UsePlot=1; 
UseSave=1; 

if(UsePlot==1)  
    %% plot variables 
   figure(); 
    fig = plot(data(:,1),data(:,2),data(:,1),data(:,3),'Linewidth',2);
    xlabel('Wavelength [nm] ','Fontsize',lfsize,'fontname',fonm);
    ylabel('Permittivity \epsilon ','Fontsize',lfsize,'fontname',fonm);
    %%%ylim([0 100]);
    %%%title(text(sprintf('%s',filename)),'fontsize',lfsize,'fontname',fonm);  
    %%%savefigname = strcat('filename',extension);
    %%%saveas(fig,savefigname);
end

data_smooth = zeros(size(data));
data_smooth(:,1) = data(:,1); %%data_smooth(:,1) is wavelength 
data_smooth(:,3) = data(:,3); 
%    data(:,1) = 2*pi*1000./data(:,1);
data(:,1) = 2*pi*10^9*c_light./data(:,1);  %%data(:,1) is omega 
data(:,2) = smooth(data(:,2),0.1,'loess'); 
data_smooth(:,2) = data(:,2); 

if(UseSave==1)
    fid = fopen(outname,'w');
    formatSpec= '%8.6f    %8.6f+%8.6fi\n';
%    fprintf(fid,formatSpec,data);
    for(ii=1:length(data))
        fprintf(fid,formatSpec, data(ii,:));
        disp(data(ii,:));
    end
    fclose(fid);
    
fid = fopen('wvnm_eps_raw/Au_Palik_wvnm_eps_smooth.dat','w');
    formatSpec= '%8.6f    %8.6f+%8.6fi\n';
    for(ii=1:length(data_smooth))
        fprintf(fid,formatSpec, data_smooth(ii,:));
        disp(data(ii,:));
    end
    fclose(fid);
    
end

if(UsePlot==1)  
    %% plot variables 
    fonm = 'arial';
    lfsize = 15;
    sfsize = 12;
    %%colorlist= 'rgbkmcrgbkmc';
    extension = '.jpg';
    fig =figure(); 
     plot(data(:,1),data(:,2),data(:,1),data(:,3),'Linewidth',2);
    xlabel('Angular frequency [\omega] ','Fontsize',lfsize,'fontname',fonm);
    ylabel('Permittivity \epsilon ','Fontsize',lfsize,'fontname',fonm);
    %%ylim([0 100]);
    %title(text(sprintf('%s',filename)),'fontsize',lfsize,'fontname',fonm);  
    savefigname = strcat(filename,extension);
    saveas(fig,savefigname);
end


cd (basedir);