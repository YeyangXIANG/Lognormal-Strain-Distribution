%% Strain statistical analysis

%%  Clear and close all

clear all; close all; clc
warning off

%% Load Nccor data

% Read file
[DataName, DataPath] = uigetfile('*.mat', 'Choose Ncorr data file');
cd(DataPath);
load([DataPath DataName]);

% Load DIC data from Ncorr

for i = 1:length(data_dic_save.strains)
    exx{i} = data_dic_save.strains(i).plot_exx_ref_formatted;
    exy{i} = data_dic_save.strains(i).plot_exy_ref_formatted;
    eyy{i} = data_dic_save.strains(i).plot_eyy_ref_formatted;
end

% Choose the strain type, exx, exy or eyy
prompt = ['Choose, 1: Strain exx,    2: Strain exy,    3: Strain eyy'];
titlename = ['Input the number'];
answer = inputdlg(prompt,titlename);
if answer{1} == '1'
    StrainType = 'exx';
elseif answer{1} == '2'
    StrainType = 'exy';
elseif answer{1} == '3'
    StrainType = 'eyy';
end

try
    disp(['The strain ' StrainType ' will be analyzed'])
catch
    disp('Please enter the right number')
end
Strain = eval(StrainType);

%% Read image files

if isempty(dir('*.JPG'))
    ImageFiles = dir('*.tif');
else
    ImageFiles = dir('*.JPG');
end

I_tmp = imread(ImageFiles(1).name);

%% Create folder to save results

disp(['Folder ' StrainType ' is created to save results'])
mkdir(StrainType);
cd([DataPath '\' StrainType])

%% Delete region should not be analyzed (only for Ncorr data)

% Choose region
figure(99)
plotmap(Strain{length(Strain)});
title('Click four corners in the figure')
trim_image = ginput(4);
r1 = mink(trim_image,2); r1 = r1(2,:);
r2 = maxk(trim_image,2); r2 = r2(2,:);
close(99)

for i = 1:length(Strain)
    Strain{i} = Strain{i}(r1(2):r2(2),r1(1):r2(1));
end

% Replace the bad points
for i = 1:length(Strain)
    ind = find(Strain{i} == 0);
    Strain{i}(ind) = NaN;
    Strain{i} = fillmissing(Strain{i});
end

%% Plot strain maps (only for Ncorr data)

for i = 1:length(Strain)
    figure(100)
    plotmap(Strain{i});
    title(['Strain ' StrainType ', image ' ImageFiles(i+1).name(1:end-4)],'interpreter','latex');
    axis off
    print(gcf,['Strain ' StrainType ', image ' ImageFiles(i+1).name(1:end-4)],'-dpng','-r400');
    saveas(gcf,['Strain ' StrainType ', image ' ImageFiles(i+1).name(1:end-4)],'fig');
    close(100)
end

%% Calculate strain, std directly from DIC data

for i = 1:length(Strain)
    strain_DIC_median(i) = median(Strain{i},'all');
    strain_DIC_mean(i) = mean(Strain{i},'all');
    strain_std(i) = std(Strain{i},[],'all');
end

%% Fit log-normal distribution

for i = 1:length(Strain)
    
    % get data > 0
    ind_0 = find(Strain{i} > 0);
    if isempty(ind_0)
        continue;
    end
    e = reshape(Strain{i}(ind_0),[],1);
    
    % divide data into 100 parts
    figure(109)
    h=histogram(e,100,'Normalization','pdf');
    x = h.BinEdges;
    x = x(2:end);
    y = h.Values;
    
    % fit log-normal
    pd(i) = fitdist(e,'Lognormal');
    
    % confidence intervals with significance level 95%
    ci{i} = paramci(pd(i));
    
    % R square
    y_fit = pdf(pd(i),x);
    R2(i) = CalR2(y,y_fit);
    
    % Sigma and median of strain
    sigma(i) = pd(i).sigma;
%     strain_median(i) = exp(pd(i).mu);
    strain_median(i) = exp(pd(i).mu);
    
    % plot figures
    %set(h,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5]);
    hold on;
    plot(x,y_fit,'r-','linewidth',2);
    hold off
    title('Fit result of lognormal distribution')
    xlabel('Strain')
    ylabel('Probability density')
    text(0.65,0.4,['\sigma = ' num2str(sigma(i))],'color','b','fontsize',15,'units','normalized');
    text(0.65,0.5,['\epsilon_c = ' num2str(strain_median(i))],'color','b','fontsize',15,'units','normalized');
    text(0.65,0.6,['R^2 = ' num2str(R2(i))],'color','b','fontsize',15,'units','normalized');
    % xlim([min(x) max(x)]);
    boxaxes;
% Nccorr
% 	print(gcf,['Lognormal Strain ' StrainType ' image ' ImageFiles(i).name(1:end-4)],'-dpng','-r400');
%     saveas(gcf,['Lognormal Strain ' StrainType ' image ' ImageFiles(i).name(1:end-4)],'fig');
    print(gcf,['Lognormal Strain ' StrainType ', image ' ImageFiles(i+1).name(1:end-4)],'-dpng','-r400');
    saveas(gcf,['Lognormal Strain ' StrainType ', image ' ImageFiles(i+1).name(1:end-4)],'fig');
    close(109)
end

%%%% Calculate error bar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(ci)
    try
        neg_sigma(i) = sigma(i) - ci{i}(1,2);
        pos_sigma(i) = ci{i}(2,2) - sigma(i);
        neg_strain(i) = strain_median(i) - exp(ci{i}(1,1));
        pos_strain(i) = exp(ci{i}(2,1)) - strain_median(i);
    catch
        continue
    end
end

display('Finish fit lognormal distribution')

%% Folder to save result's summary

mkdir('Results');
cd([DataPath '\' StrainType '\Results'])

%% Plot evolution of results

%%%% Define strain of one image %%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = strain_DIC_mean;

%%%% sigma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% errorbar(x,sigma,neg_sigma,pos_sigma,'o')
plot(x,sigma,'o','linewidth',1.5)
ylim([0,max(sigma)*1.2])
xlabel('Mean strain based on DIC')
ylabel('\sigma of lognormal distribution')
ax.FontName = 'Times New Roman';
boxaxes;
print(gcf,['Results sigma ' StrainType],'-dpng','-r400');
saveas(gcf,['Results sigma ' StrainType],'fig');

%%%% Strain based on lognormal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(1:length(strain_median),strain_median,'o','linewidth',1.5);
xlabel('Image number')
ylabel('Strain')
title('Strain based on lognormal distribution')
ax.FontName = 'Times New Roman';
boxaxes;
print(gcf,['Strain ' StrainType '_Lognormal'],'-dpng','-r400');
saveas(gcf,['Strain ' StrainType '_Lognormal'],'fig');

%%%% Strain based on DIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(1:length(strain_DIC_median),strain_DIC_median,'ro','linewidth',1.5);
xlabel('Image number')
ylabel('Strain')
title('Strain based on DIC data')
ax.FontName = 'Times New Roman';
boxaxes;
print(gcf,['Strain ' StrainType '_DIC'],'-dpng','-r400');
saveas(gcf,['Strain ' StrainType '_DIC'],'fig');

%%%% Compare different strain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Length_s = 1:length(strain_DIC_median);
figure
% errorbar(x,strain_median,neg_strain,pos_strain,'o')
plot(Length_s,strain_median,'o',Length_s,strain_DIC_median,'o',Length_s,strain_DIC_mean,'o','linewidth',1.5);
if strain_DIC_median(1) > 0
    ylim([0,max([strain_median;strain_DIC_median;strain_DIC_mean],[],'all')*1.2])
end
xlabel('Image number')
ylabel('Strain')
title('Compare strain')
ax.FontName = 'Times New Roman';
boxaxes;
legend('Strain_{median}^{lognormal}','Strain_{median}^{DIC}','Strain_{mean}^{DIC}','location','southeast','fontsize',13)
legend('boxoff')
print(gcf,['Strain ' StrainType '_Compare'],'-dpng','-r400');
saveas(gcf,['Strain ' StrainType '_Compare'],'fig');

%%%% std %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(x,strain_std,'o','linewidth',1.5)
ylim([0,max(strain_std)*1.2])
xlabel('Mean strain based on DIC')
ylabel('std of strain')
ax.FontName = 'Times New Roman';
boxaxes;
print(gcf,['Results std ' StrainType],'-dpng','-r400');
saveas(gcf,['Results std ' StrainType],'fig');

display('Finish plot results')

%% Save results

SaveName = ['Strain_Statistics_' StrainType '.mat'];
save(SaveName,'Strain','StrainType','strain_median','sigma','R2',...
    'strain_DIC_median','strain_DIC_mean','strain_std')

warning on

display('Finish all')

%% Local functions
% Calculate delta_r along x direction

function dx = ddx(E)
Nx=size(E,2);
Nx=fix(Nx/2);
dx1=(1:Nx-1);dx2=2.^(1:0.25:log2(length(dx1)-1));dx3=round(dx2);sub=[];
if dx3(length(dx3))<Nx-1
    dx4=[1 dx3 Nx-1];
else
    dx4=[1 dx3];
end
for i=1:length(dx4)-1
    if dx4(i)==dx4(i+1)
        sub=[sub i];
    end
end
dx=dx4(setdiff((1:length(dx4)),sub));
end

%% Local function

function [deltax,deltah_square,stddeltaz]=hurst_dh_bandwidth3(dx,z)

[nrow,ncolumn]=size(z);

deltax=[];
deltah_square=[];

stddeltaz=[];
N=[];
for inddxcurrent=1:length(dx)
   dxcurrent=dx(inddxcurrent);
   zcurrent1=z(:,1:end-dxcurrent);
   zcurrent2=z(:,1+dxcurrent:end);
   dzcurrent2=(zcurrent2-zcurrent1).^2;
   deltax=[deltax;dxcurrent];
   deltah_square=[deltah_square;mean2((dzcurrent2))];
   stddeltaz=[stddeltaz;std2(dzcurrent2)];
   N=[N;prod(size(dzcurrent2))];
end

stddeltaz=stddeltaz./sqrt(N);
end

%% Local function

function R2 = CalR2(z,z_est)
% calcuateR2 Cacluate R-squared
% R2 = calcuateR2(z,z_est) takes two inputs - The observed data x and its
% estimation z_est (may be from a regression or other model), and then
% compute the R-squared value a measure of goodness of fit. R2 = 0
% corresponds to the worst fit, whereas R2 = 1 corresponds to the best fit.

r = z-z_est;
normr = norm(r);
SSE = normr.^2;
SST = norm(z-mean(z))^2;
R2 = 1 - SSE/SST;
end

%% Local function

function [a,b]=boxaxes
% get handle to current axes
% a = gca;
% % % get the original xlim, ylim
% % xlim_ori = get(a,'XLim');
% % ylim_ori = get(a,'YLim');
% % set box property to off and remove background color
% set(a,'box','off','color','none')
% % create new, empty axes with box but without ticks
% b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'xlim',get(a,'XLim'),'ylim',get(a,'YLim'));
% % set original axes as active
% axes(a)
% % link axes in case of zooming
% linkaxes([a b])

box off;
a = gca;
x = get(a,'XLim');
y = get(a,'YLim');
line([x(1) x(2)],[y(2) y(2)],'Color','k','linewidth',0.54)
line([x(2) x(2)],[y(1) y(2)],'Color','k','linewidth',0.54)
xlim(x)
ylim(y)
end

%% Local function

function plotmap(M,transparency)

carte_a_tracer=M;
Z=reshape(carte_a_tracer,1,size(carte_a_tracer,1) ...
   *size(carte_a_tracer,2));
stdZ=std(Z);
meanZ=mean(Z);
Zmin=max([min(Z) meanZ-4*stdZ]);
Zmax=min([max(Z) meanZ+4*stdZ]);
I=find(Z<=Zmin);Z(I)=Zmin.*ones(1,length(I));
I=find(Z>=Zmax);Z(I)=Zmax*ones(1,length(I));
Z=reshape(Z,size(carte_a_tracer,1),size(carte_a_tracer,2));

if (nargin<2)
    transparency = 1;
end

colormap(jet(256));
image(Z,'CDataMapping','scaled','AlphaData',transparency);
colorbar
end

%% Local function

function B=fillmissing(A)

[n,m]=size(A);
A=A(:);
nm=n*m;
k=isnan(A(:));

% list the nodes which are known, and which will
% be interpolated
nan_list=find(k);
known_list=find(~k);

% how many nans overall
nan_count=length(nan_list);

% convert NaN indices to (r,c) form
% nan_list==find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc]=ind2sub([n,m],nan_list);

% both forms of index in one array:
% column 1 == unrolled index
% column 2 == row index
% column 3 == column index
nan_list=[nan_list,nr,nc];


   hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
  hv_springs=[];
  for i=1:4
    hvs=nan_list+repmat(hv_list(i,:),nan_count,1);
    k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
    hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
  end

  % delete replicate springs
  hv_springs=unique(sort(hv_springs,2),'rows');
  
  % build sparse matrix of connections, springs
  % connecting diagonal neighbors are weaker than
  % the horizontal and vertical springs
  nhv=size(hv_springs,1);
  springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
     repmat([1 -1],nhv,1),nhv,nm);
  
  % eliminate knowns
  rhs=-springs(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;

B=reshape(B,n,m);
function neighbors_list=identify_neighbors(n,m,nan_list,talks_to)

if ~isempty(nan_list)
  % use the definition of a neighbor in talks_to
  nan_count=size(nan_list,1);
  talk_count=size(talks_to,1);
  
  nn=zeros(nan_count*talk_count,2);
  j=[1,nan_count];
  for i=1:talk_count
    nn(j(1):j(2),:)=nan_list(:,2:3) + ...
        repmat(talks_to(i,:),nan_count,1);
    j=j+nan_count;
  end
  
  % drop those nodes which fall outside the bounds of the
  % original array
  L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m); 
  nn(L,:)=[];
  
  % form the same format 3 column array as nan_list
  neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];
  
  % delete replicates in the neighbors list
  neighbors_list=unique(neighbors_list,'rows');
  
  % and delete those which are also in the list of NaNs.
  neighbors_list=setdiff(neighbors_list,nan_list,'rows');
  
else
  neighbors_list=[];
end


end
end