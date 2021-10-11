 % First order and total effect indices for a given
% model computed with Extended Fourier Amplitude
% Sensitivity Test (EFAST).
% Andrea Saltelli, Stefano Tarantola and Karen Chan.
% 1999. % "A quantitative model-independent method for global
% sensitivity analysis of model output". % Technometrics 41:3956
clc
clear;
close all;
tic
global y_var_label
poolobj = gcp;
addAttachedFiles(poolobj,{'Parameter_settings_EFAST.m',...
    'CVmethod.m','efast_sd.m','efast_ttest.m',...
    'parameterdist.m','INS1Epathway.m','SETFREQ.m','Run_simulation_sensitivity.m','simulation_protocol_sensitivity.m'}); % 'simulation_protocol_fitting_Spegel2015.m'})
runName = ['20210601_GrahamEfast'];
mkdir(runName);

%% PARAMETERS AND ODE SETTINGS (they are included in the following file)
Parameter_settings_EFAST;
%% INPUT
NR = 5; %: no. of search curves - RESAMPLING
k = length(efast_var); % # of input factors (parameters varied) + dummy parameter
NS = 125; % # of samples per search curve
wantedN=NS*k*NR; % wanted no. of sample points
% OUTPUT
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of k frequencies
% OMi : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output
MI = 4; %: maximum number of fourier coefficients
% that may be retained in calculating the partial
% variances without interferences between the
% assigned frequencies
% Computation of the frequency for the group
% of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
        '65 per factor.\n']);
    return;
end
%diary my_data.txt
%% Pre-allocation of the output matrix Y
%% Y will save only the points of interest specified in
%% the vector time_points
Y(NS,1,length(y_var_label),length(pmin),NR)=0;  % pre-allocation

% Loop over k parameters (input factors)
for i=1:k % i=# of replications (or blocks)
    % Algorithm for selecting the set of frequencies.
    % OMci(i), i=1:k-1, contains the set of frequencies
    % to be used by the complementary group.
    OMci = SETFREQ(k,OMi/2/MI,i);
    % Loop over the NR search curves.
    for L=1:NR
        % Setting the vector of frequencies OM
        % for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar
        % variable S and the coordinates
        % {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1
        % Transform distributions from standard uniform to general.
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmin,pmax,[],[],NS,'unif'); %%this is what assigns 'our' values rather than 0:1 dist
        % Do the NS model evaluations.
        %paramsdidnotwork = zeros(NS,length(pmin));
        parfor run_num=1:NS
            % This gives the output for
            % keeps track of [parameter run NR]
            [i run_num L]
            % ODE system file
            % @pathway_ode.m
            % ODE solver call t,y,params,X,run_num
             try
              % [t,y] = Run_simulation_sensitivity(X(:,:,i,L),run_num);
               [t,y] = Run_simulation_sensitivity(X(:,:,i,L),run_num);

             catch
                err = 1;
                fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
                y = zeros(1,112);
             end
           Y(run_num,:,:,i,L)= y(:,:);% time points of interest for the US analysis
        end %run_num=1:NS
        %dfile = ['for_k_' num2str(i) '_' num2str(L) '.mat'];
        %save(dfile,'paramsdidnotwork','-v7.3')
    end % L=1:NR
end % i=1:k
%diary off
save([runName '/Model_efast.mat']);

%CALCULATE Si AND STi for each resample (1,2,...,NR) [ranges]
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,1,1:length(y_var_label));
save([runName '/Si.mat'],'Si')
save([runName '/Sti.mat'],'Sti')
save([runName '/rangeSi.mat'],'rangeSi')
save([runName '/rangeSti.mat'],'rangeSti')

for SensCalc = 1:length(y_var_label)
    ind = SensCalc;
    [CVsi(:,:,SensCalc)  CVsti(:,:,SensCalc)]=CVmethod(Si, rangeSi,Sti,rangeSti,ind);
    save([runName '/CVsi.mat'],'CVsi')
    save([runName '/CVsti.mat'],'CVsti')
    %T-test on Si and STi for Viral load (variable 4)
    %[eFASToutputs{SensCalc}, Si_out(:,:,SensCalc), p_Si_out(:,:,SensCalc), Sti_out(:,:,SensCalc), p_Sti_out(:,:,SensCalc)] = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),efast_var,ind,y_var_label,0.05);
    [eFASToutputs{SensCalc}, Si_out(:,:,SensCalc), p_Si_out(:,:,SensCalc), Sti_out(:,:,SensCalc), p_Sti_out(:,:,SensCalc)] = efast_ttest(Si,rangeSi,Sti,rangeSti,1,efast_var,ind,y_var_label,0.05);    
    save([runName '/eFASToutputs.mat'],'eFASToutputs')
    save([runName '/Si_out_data.mat'],'Si_out')
    save([runName '/p_Si_out_data.mat'],'p_Si_out')
    save([runName '/Sti_out_data.mat'],'Sti_out')
    save([runName '/p_Sti_out_data.mat'],'p_Sti_out')
end


p_Si_out = squeeze(p_Si_out)';
p_Sti_out = squeeze(p_Sti_out)';
Si_out = squeeze(Si_out)';
Sti_out = squeeze(Sti_out)';

dummy_Si_out = Si_out(:,end);
dummy_Sti_out = Sti_out(:,end);

Si_out = Si_out(:,1:end-1);
Sti_out = Sti_out(:,1:end-1);


[a b] = size(Si_out);

Si_out(find(p_Si_out>0.05)) = NaN;
Sti_out(find(p_Sti_out>0.05)) = NaN;

for i = 1:a
    Si_out(i,find(Si_out(i,:)<=dummy_Si_out(i))) = NaN;
    Sti_out(i,find(Sti_out(i,:)<=dummy_Sti_out(i))) = NaN;
end

%Labeling
outputs =  {'Graham_3PG_2_1'
'Graham_ADP_2_1'
'Graham_AKG_2_1'
'Graham_AMP_2_1'
'Graham_ATP_2_1'
'Graham_CitICit_2_1'
'Graham_FUM_2_1'
'Graham_G3P_2_1'
'Graham_GDP_2_1'
'Graham_GLC_2_1'
'Graham_GSH_2_1'
'Graham_GSSG_2_1'
'Graham_G6P_2_1'
'Graham_ALA_2_1'
'Graham_GLUT_2_1'
'Graham_GLN_2_1'
'Graham_LAC_2_1'
'Graham_MAL_2_1'
'Graham_NADP_2_1'
'Graham_NADPH_2_1'
'Graham_Ru5P_2_1'
'Graham_SUC_2_1'
'Graham_3PG_16_1'
'Graham_ADP_16_1'
'Graham_AKG_16_1'
'Graham_AMP_16_1'
'Graham_ATP_16_1'
'Graham_CitICit_16_1'
'Graham_FUM_16_1'
'Graham_G3P_16_1'
'Graham_GDP_16_1'
'Graham_GLC_16_1'
'Graham_GSH_16_1'
'Graham_GSSG_16_1'
'Graham_G6P_16_1'
'Graham_ALA_16_1'
'Graham_GLUT_16_1'
'Graham_GLN_16_1'
'Graham_LAC_16_1'
'Graham_MAL_16_1'
'Graham_NADP_16_1'
'Graham_NADPH_16_1'
'Graham_Ru5P_16_1'
'Graham_SUC_16_1'
'Graham_3PG_25_1'
'Graham_ADP_25_1'
'Graham_AKG_25_1'
'Graham_AMP_25_1'
'Graham_ATP_25_1'
'Graham_CitICit_25_1'
'Graham_FUM_25_1'
'Graham_G3P_25_1'
'Graham_GDP_25_1'
'Graham_GLC_25_1'
'Graham_GSH_25_1'
'Graham_GSSG_25_1'
'Graham_G6P_25_1'
'Graham_ALA_25_1'
'Graham_GLUT_25_1'
'Graham_GLN_25_1'
'Graham_LAC_25_1'
'Graham_MAL_25_1'
'Graham_NADP_25_1'
'Graham_NADPH_25_1'
'Graham_Ru5P_25_1'
'Graham_SUC_25_1'
'Graham_3PG_16_2'
'Graham_ADP_16_2'
'Graham_AKG_16_2'
'Graham_AMP_16_2'
'Graham_ATP_16_2'
'Graham_CitICit_16_2'
'Graham_FUM_16_2'
'Graham_G3P_16_2'
'Graham_GDP_16_2'
'Graham_GLC_16_2'
'Graham_GSH_16_2'
'Graham_GSSG_16_2'
'Graham_G6P_16_2'
'Graham_ALA_16_2'
'Graham_GLUT_16_2'
'Graham_GLN_16_2'
'Graham_LAC_16_2'
'Graham_MAL_16_2'
'Graham_NADP_16_2'
'Graham_NADPH_16_2'
'Graham_Ru5P_16_2'
'Graham_SUC_16_2'
'Graham_3PG_25_2'
'Graham_ADP_25_2'
'Graham_AKG_25_2'
'Graham_AMP_25_2'
'Graham_ATP_25_2'
'Graham_CitICit_25_2'
'Graham_FUM_25_2'
'Graham_G3P_25_2'
'Graham_GDP_25_2'
'Graham_GLC_25_2'
'Graham_GSH_25_2'
'Graham_GSSG_25_2'
'Graham_G6P_25_2'
'Graham_ALA_25_2'
'Graham_GLUT_25_2'
'Graham_GLN_25_2'
'Graham_LAC_25_2'
'Graham_MAL_25_2'
'Graham_NADP_25_2'
'Graham_NADPH_25_2'
'Graham_Ru5P_25_2'
'Graham_SUC_25_2'
'Graham_ASP_16_2'
'Graham_ASP_25_2'};
% 
% outputs={'Spegel_3 2pg'
% 'Spegel_3 3pg'
% 'Spegel_3 akg'
% 'Spegel_3 ala'
% 'Spegel_3 asp'
% 'Spegel_3 cit'
% 'Spegel_3 fum'
% 'Spegel_3 g3p'
% 'Spegel_3 lac'
% 'Spegel_3 mal'
% 'Spegel_3 pep'
% 'Spegel_3 pyr'
% 'Spegel_3 r5p'
% 'Spegel_3 suc'
% 'Spegel_6 2pg'
% 'Spegel_6 3pg'
% 'Spegel_6 akg'
% 'Spegel_6 ala'
% 'Spegel_6 asp'
% 'Spegel_6 cit'
% 'Spegel_6 fum'
% 'Spegel_6 g3p'
% 'Spegel_6 lac'
% 'Spegel_6 mal'
% 'Spegel_6 pep'
% 'Spegel_6 pyr'
% 'Spegel_6 r5p'
% 'Spegel_6 suc'
% 'Spegel_10 2pg'
% 'Spegel_10 3pg'
% 'Spegel_10 akg'
% 'Spegel_10 ala'
% 'Spegel_10 asp'
% 'Spegel_10 cit'
% 'Spegel_10 fum'
% 'Spegel_10 g3p'
% 'Spegel_10 lac'
% 'Spegel_10 mal'
% 'Spegel_10 pep'
% 'Spegel_10 pyr'
% 'Spegel_10 r5p'
% 'Spegel_10 suc'
% 'Spegel_15 2pg'
% 'Spegel_15 3pg'
% 'Spegel_15 akg'
% 'Spegel_15 ala'
% 'Spegel_15 asp'
% 'Spegel_15 cit'
% 'Spegel_15 fum'
% 'Spegel_15 g3p'
% 'Spegel_15 lac'
% 'Spegel_15 mal'
% 'Spegel_15 pep'
% 'Spegel_15 pyr'
% 'Spegel_15 r5p'
% 'Spegel_15 suc '
% 'Malmgren AKG '
% 'Malmgren ALA '
% 'Malmgren ASP '
% 'Malmgren CIT '
% 'Malmgren FUM '
% 'Malmgren G6P '
% 'Malmgren GLC '
% 'Malmgren glut'
% 'Malmgren G3P '
% 'Malmgren ICIT'
% 'Malmgren LAC '
% 'Malmgren MAL '
% 'Malmgren PYR '
% 'Malmgren SUC'
% 'Goehring FRU'};


labels_inputs = {'Vf_ glut'
    'Vf_ gk'
    'Vf_ hpi '
    'Vr_ hpi '
    'Vf_ pfk1'
    'Vf_ aldo'
    'Vr_ aldo'
    'Vf_ tpi '
    'Vf_ gapdh '
    'Vf_ pgk '
    'Vf_ pgam'
    'Vf_ eno '
    'Vf_ pyk '
    'Vf_ ldh '
    'Vr_ ldh '
    'Vf_ ak'
    'Vf_ atpase'
    'Vr_ atpase'
    'Vf_ ox'
    'Vf_ mct1'
    'Vf_ g6pd'
    'Vr_ g6pd'
    'Vf_ 6pgdh '
    'Vr_ 6pgdh '
    'Vf_ rpe '
    'Vf_ rpi '
    'Vf_ prpps '
    'Vf_ tk1 '
    'Vr_ tk1 '
    'Vf_ tk2 '
    'Vf_ ta'
    'Vf_ gpx '
    'Vf_ gssgr '
    'Vr_ gssgr '
    'Vf_ pdh '
    'Vf_ cs'
    'Vf_ acon'
    'Vf_ idh '
    'Vf_ akgd'
    'Vf_ s '
    'Vf_ sdh '
    'Vf_ fum '
    'Vf_ mdh2'
    'Vf_ got2'
    'Vf_ mdh1'
    'Vf_ got1'
    'Vf_ akgmal'
    'Vf_ aspglu'
    'Vf_ pyrh'
    'Vf_ citmal'
    'Vf_ malpi '
    'Vf_ gluh'
    'Vf_ cly '
    'Vf_ malic '
    'Vf_ cmalic'
    'Vf_ pc'
    'Vf_ gls '
    'Vf_ gdh '
    'Vf_ gpt '
    'Vf_ asct2 '
    'f1_ asct2 '
    'Kgluout_ asct2'
    'Kgluin_ asct2 '
    'f2_ asct2 '
    'Keq1_ asct2 '
    'Kgluout1_ asct2 '
    'Kgluin1_ asct2'
    'Vf_ aconitase '
    'Vf_ cIDH'
    'Vf_ isocitmal '
    'Glu_out'
    'gk_ K1GLC'
    'gk_ K1ATP'
    'E_ aldr'
    'Vm_ SoDH '
    'kfruT '
    'Ct_ PyP'
    'Ct_ Pyr'
    'k1_ aldr'
    'Vm_r_SoDH'
    'Vf_hk'
    'dummy'}';
figure()
ny = length(outputs);
nx = length(labels_inputs);
%padded = padarray((Sti_out),[1,1],'post');
%save([runName '/Padded.mat'],'padded')
%p = pcolor(padded);
p = pcolor(Sti_out);

%p = pcolor(meanPad);

colorDepth = 1000;
% colormap(newColorMap(colorDepth));
%colormap default
colormap(flipud(bone))
%set(p, 'Edgecolor','black','LineWidth',2)
set(gca,'xgrid', 'off', 'ygrid', 'off', 'gridlinestyle', '-', 'Xcolor','k', 'Ycolor', 'k','LineWidth',3);
set(gca,'FontSize',20, 'FontWeight','bold','Fontname','Arial')
pbaspect([nx ny 1])
% axis square;
% axis tight
colorbar
caxis ([0  1]);
title([runName 'Efast Run '])
set(gca,'XTick',1:1:b);
set(gca,'YTick',1:1:a);
set(gcf,'color','white')
ax = gca;
XTick = get(ax, 'XTick');
XTickLabel = get(ax, 'XTickLabel');
set(ax,'XTick',XTick+0.5)
set(ax,'XTickLabel',XTickLabel)
YTick = get(ax, 'YTick');
YTickLabel = get(ax, 'YTickLabel');
set(gca,'TickLength',[ 0 0 ])
set(ax,'YTick',YTick+0.5)
set(ax,'YTickLabel',YTickLabel)
set(gca,'XTickLabel',labels_inputs)
set(gca,'YTickLabel',outputs)
set(gca,'XTickLabelRotation',90)
set(gcf, 'renderer', 'painters');
set(gcf, 'color', 'white');
f = get(gca,'title');
set(f,'FontSize',20,'FontWeight','bold')
filename = [runName '/Sti.fig'];
savefig(filename);
toc
delete(gcp)
