function [timepoints, observables_out] = Run_simulation_fitting(pos)
%params = zeros(59,1);
params = pos;

%% Solving the differntial equations
[timepoints, observables_out] = simulation_protocol_fitting(params);



%% Glycolysis
% params(1,1)         = pos(1);       %Vf_glut = 0.028
% params(2,1)         = pos(2);       %Vf_hk = 0.041
% params(3,1)         = pos(3);       %Vf_hpi = 0.28
% params(4,1)         = pos(4);       %Vf_pfk1 = 0.022
% params(5,1)         = pos(5);       %Vf_aldo = 0.08
% params(6,1)         = pos(6);       %Vf_tpi = 3.4              
% params(7,1)         = pos(7);       %Vf_gapdh = 0.331
% params(8,1)         = pos(8);       %Vr_gapdh = 0.413
% params(9,1)         = pos(9);       %Vf_pgk = 8.7
% params(10,1)        = pos(10);      %Vf_pgam = 0.94
% params(11,1)        = pos(11);      %Vf_eno = 0.34
% params(12,1)        = pos(12);      %Vf_pyk = 0.087
% params(13,1)        = pos(13);      %Vf_ldh = 0.468
% params(14,1)        = pos(14);      %Vr_ldh = 0.074
% params(15,1)        = pos(15);      %Vf_ak = 200/6
% params(16,1)        = pos(16);      %Vf_atpase = 39/6
% params(17,1)        = pos(17);      %Vr_atpase = 2.1179e-14
% params(18,1)        = pos(18);      %Vf_ox = 0.140
% params(19,1)        = pos(19);      %Vf_mct1 = 0.03
% %% PPP pathway
% params(20,1)        = pos(20);      %Vf_g6pd = (4.39e+20)/60
% params(21,1)        = pos(21);      %Vf_6pgdh = (1.83e+19)/60
% params(22,1)        = pos(22);      %Vf_rpe = (4.642e+04)/60
% params(23,1)        = pos(23);      %Vf_rpi = (2.56e+04)/60
% params(24,1)        = pos(24);      %Vf_prpps = 25.3/60
% params(25,1)        = pos(25);      %Vf_tk1 = (5.32e+18)/60
% params(26,1)        = pos(26);      %Vr_tk1 = (2.58e+18)/60
% params(27,1)        = pos(27);      %Vf_tk2 = (3.64e+19)/60
% params(28,1)        = pos(28);      %Vr_tk2 = (1.23e+18)/60
% params(29,1)        = pos(29);      %Vf_ta = (1.23e+18)/60 
% params(30,1)        = pos(30);      %Vr_ta = (7.53e+19)/60 
% params(31,1)        = pos(31);      %Vf_gpx = 1.56e+04/60
% params(32,1)        = pos(32);      %Vf_gssgr = 5.5e+044/60
% params(33,1)        = pos(33);      %Vr_gssgr = 1.05e+40/60
% %% TCA cycle reactions
% params(34,1)        = pos(34);      %Vf_pdh = 189.7/60
% params(35,1)        = pos(35);      %Vf_cs = 1804/60
% params(36,1)        = pos(36);      %Vf_acon = 100000/60
% params(37,1)        = pos(37);      %Vf_idh = 27500/60
% params(38,1)        = pos(38);      %Vf_akgd = 2495/60
% params(39,1)        = pos(39);      %Vf_s = 362/60
% params(40,1)        = pos(40);      %Vf_sdh = 5.81e+04/60
% params(41,1)        = pos(41);      %Vf_fum = 2.21e+04/60
% params(42,1)        = pos(42);      %Vf_mdh2 = 3.53e+08/60
% params(43,1)        = pos(43);      %Vf_got2 = 3.87e+06/60 
% params(44,1)        = pos(44);      %Vf_mdh1 = 3.59e+07/60
% params(45,1)        = pos(45);      %Vf_got1 = 2.36e+03/60
% params(46,1)        = pos(46);      %Vmf_akgmal = 3.19e+06/60
% params(47,1)        = pos(47);      %Vmf_aspglu = 2.49e+04/60
% params(48,1)        = pos(48);      %Vm_pyrh = 10e+13/60
% params(49,1)        = pos(49);      %Vm_citmal = 296.6/60
% params(50,1)        = pos(50);      %Vm_malpi = 17.3/60
% params(51,1)        = pos(51);      %Vm_gluh = 3.87e+08/60
% params(52,1)        = pos(52);      %Vm_cly = 17.5/60
% params(53,1)        = pos(53);      %Vm_malic = 46.3/60
% params(54,1)        = pos(54);      %Vf_cmalic   = 174/60
% params(55,1)        = pos(55);      %V_pc = 718.2/60
% params(56,1)        = pos(56);      %Vm_gls = 38.8/60
% params(57,1)        = pos(57);      %Vf_gdh = 5550/60
% params(58,1)        = pos(58);      %Vf_gpt = 66.4/60
% params(59,1)        = pos(59);      %Vmax_asct2 = 0.02
 