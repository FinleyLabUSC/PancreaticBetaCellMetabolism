
%% All Predictions
function [observables] = AllPredictions_CollectedTimeCourse(p)
    params = p;
    
    %A = [0.616000000000000;3;0.726154480000000;0.108367000000000;0.659000000000000;0.00700000000000000;18.5478285800000;0.00191000000000000;1;0.299000000000000;1.00064013300000;2.36164831900000;0.00211000000000000;0.00423000000000000;77.3800000000000;0.00261149000000000;0.00963821300000000;0.00344545400000000;0.00143820600000000;0.0205015600000000;0.0759752810000000;0.00122000000000000;0.170000000000000;6.29000000000000e-06;0.102500000000000;0.0630000000000000;0.360000000000000;0.420000000000000;0.400000000000000;0.294100000000000;2.95000000000000;0.0650000000000000;0.500000000000000;0.00450000000000000;1.60000000000000;3.03000000000000;4.65000000000000e-05;2.68000000000000;0.00450000000000000;1.97000000000000e-05;0.225000000000000;4.72000000000000e-05;4.44750337500000;0.0720000000000000;0.0720000000000000;8.44445922300000;0.0240000000000000;0.450000000000000;4.50000000000000;0.700000000000000;4.77000000000000e-06;0.0146347570000000;1.00000000000000e-06;0.0164989410000000;0.379987684000000;0.272335568000000];
    A = [0.542125648135445;3.76168895009916;0.454448540472248;0.706990028881560;0.0156565460978362;0.478582126811771;4.39031164654454;0.168038706461350;4.21968660787895;0.0560506065370617;0.397377606983681;0.583587181302951;0.716708019845057;0.000154147902632731;0.136501997249423;0.132939511019276;0.00456400046649098;0.000733943694403398;0.00134582377890930;0.00483413989498314;0.0379963772949342;0.00192408120909374;4.23518730127469;2.45107434802339e-06;0.000771217660846397;2.09292143332501;0.288669565535384;0.00779378981985000;0.224895886981259;12.6827767477872;6.03285037234089e-05;6.77882129684873e-05;0.000309121127035145;1.74318833008063e-06;1.45680249381982e-05;2.94658145921870;3.12055381157148e-05;0.797796208874077;2.65472943228149e-07;2.85338922615892e-05;0.0240841285605586;1.05706088817784e-05;3.31601148362778;0.0170013921210288;0.00150069872532053;2.47590966892564;0.0300000000000000;2.15537701634457e-07;4.94999978446232;3.43977996188819;1.51575084594245e-06;0.0400301566013678;1.72526016292169;0.0300000000000000;1.03960000000000;0.0480000000000000];
    %A = params(73:126);
    options = odeset('RelTol',1e-12,...
        'AbsTol',1e-19,...
        'NormControl','on',...
        'Stats','off',...
        'BDF','off',...
        'NonNegative',1:56,...
        'MaxOrder', 5);
    
    %% Spegel 2013
    tspan = [0:120];
    % Afirst = A;
    initvalue = A;
    params(71,1) = 2.8;
    [~, BASE_Spegel2013] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
    try
        tspan = [0:15];
        initvalue = (BASE_Spegel2013(end,:))';
        params(71,1) = 16.7;
        [time, END_Spegel2013] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
        
        Spegel_3_2pg = END_Spegel2013(3,12) / BASE_Spegel2013(end, 12);
        Spegel_6_2pg = END_Spegel2013(6,12) / BASE_Spegel2013(end, 12);
        Spegel_10_2pg = END_Spegel2013(10,12) / BASE_Spegel2013(end, 12);
        Spegel_15_2pg = END_Spegel2013(end,12) / BASE_Spegel2013(end, 12);
        
        Spegel_3_3pg = END_Spegel2013(3,11) / BASE_Spegel2013(end, 11);
        Spegel_6_3pg = END_Spegel2013(6,11) / BASE_Spegel2013(end, 11);
        Spegel_10_3pg = END_Spegel2013(10,11) / BASE_Spegel2013(end, 11);
        Spegel_15_3pg = END_Spegel2013(end,11) / BASE_Spegel2013(end, 11);
        
        Spegel_3_akg = (END_Spegel2013(3,41)+END_Spegel2013(3,29)) / (BASE_Spegel2013(end, 41)+BASE_Spegel2013(end, 29));
        Spegel_6_akg = (END_Spegel2013(6,41)+END_Spegel2013(6,29)) / (BASE_Spegel2013(end, 41)+BASE_Spegel2013(end, 29));
        Spegel_10_akg = (END_Spegel2013(10,41)+END_Spegel2013(10,29)) / (BASE_Spegel2013(end, 41)+BASE_Spegel2013(end, 29));
        Spegel_15_akg = (END_Spegel2013(end,41)+END_Spegel2013(end,29)) / (BASE_Spegel2013(end, 41)+BASE_Spegel2013(end, 29));
        
        Spegel_3_ala = END_Spegel2013(3,50) / BASE_Spegel2013(end, 50);
        Spegel_6_ala = END_Spegel2013(6,50) / BASE_Spegel2013(end, 50);
        Spegel_10_ala = END_Spegel2013(10,50) / BASE_Spegel2013(end, 50);
        Spegel_15_ala = END_Spegel2013(end,50) / BASE_Spegel2013(end, 50);
        
        Spegel_3_asp = (END_Spegel2013(3,35)+END_Spegel2013(3,37)) / (BASE_Spegel2013(end, 35)+BASE_Spegel2013(end, 37));
        Spegel_6_asp = (END_Spegel2013(6,35)+END_Spegel2013(6,37)) / (BASE_Spegel2013(end, 35)+BASE_Spegel2013(end, 37));
        Spegel_10_asp = (END_Spegel2013(10,35)+END_Spegel2013(10,37)) / (BASE_Spegel2013(end, 35)+BASE_Spegel2013(end, 37));
        Spegel_15_asp = (END_Spegel2013(end,35)+END_Spegel2013(end,37)) / (BASE_Spegel2013(end, 35)+BASE_Spegel2013(end,37));
        
        Spegel_3_cit = (END_Spegel2013(3,27)+END_Spegel2013(3,42)) / (BASE_Spegel2013(end, 27)+BASE_Spegel2013(end, 42));
        Spegel_6_cit = (END_Spegel2013(6,27)+END_Spegel2013(6,42)) / (BASE_Spegel2013(end, 27)+BASE_Spegel2013(end, 42));
        Spegel_10_cit = (END_Spegel2013(10,27)+END_Spegel2013(10,42)) / (BASE_Spegel2013(end, 27)+BASE_Spegel2013(end, 42));
        Spegel_15_cit = (END_Spegel2013(end,27)+END_Spegel2013(end,42)) / (BASE_Spegel2013(end, 27)+BASE_Spegel2013(end,42));
        
        Spegel_3_fum = END_Spegel2013(3,32) / BASE_Spegel2013(end, 32);
        Spegel_6_fum = END_Spegel2013(6,32) / BASE_Spegel2013(end, 32);
        Spegel_10_fum = END_Spegel2013(10,32) / BASE_Spegel2013(end, 32);
        Spegel_15_fum = END_Spegel2013(end,32) / BASE_Spegel2013(end, 32);
        
        Spegel_3_g3p = END_Spegel2013(3,8) / BASE_Spegel2013(end, 8);
        Spegel_6_g3p = END_Spegel2013(6,8) / BASE_Spegel2013(end, 8);
        Spegel_10_g3p = END_Spegel2013(10,8) / BASE_Spegel2013(end, 8);
        Spegel_15_g3p = END_Spegel2013(end,8) / BASE_Spegel2013(end, 8);
        
        Spegel_3_lac = END_Spegel2013(3,15) / BASE_Spegel2013(end, 15);
        Spegel_6_lac = END_Spegel2013(6,15) / BASE_Spegel2013(end, 15);
        Spegel_10_lac = END_Spegel2013(10,15) / BASE_Spegel2013(end, 15);
        Spegel_15_lac = END_Spegel2013(end,15) / BASE_Spegel2013(end, 15);
        
        Spegel_3_mal = (END_Spegel2013(3,33)+END_Spegel2013(3,40)) / (BASE_Spegel2013(end, 33)+BASE_Spegel2013(end, 40));
        Spegel_6_mal = (END_Spegel2013(6,33)+END_Spegel2013(6,40)) / (BASE_Spegel2013(end, 33)+BASE_Spegel2013(end, 40));
        Spegel_10_mal = (END_Spegel2013(10,33)+END_Spegel2013(10,40)) / (BASE_Spegel2013(end, 33)+BASE_Spegel2013(end, 40));
        Spegel_15_mal = (END_Spegel2013(end,33)+END_Spegel2013(end,40)) / (BASE_Spegel2013(end, 33)+BASE_Spegel2013(end,40));
        
        Spegel_3_pep = END_Spegel2013(3,13) / BASE_Spegel2013(end, 13);
        Spegel_6_pep = END_Spegel2013(6,13) / BASE_Spegel2013(end, 13);
        Spegel_10_pep = END_Spegel2013(10,13) / BASE_Spegel2013(end, 13);
        Spegel_15_pep = END_Spegel2013(end,13) / BASE_Spegel2013(end, 13);
        
        Spegel_3_pyr = (END_Spegel2013(3,14)+END_Spegel2013(3,25)) / (BASE_Spegel2013(end, 14)+BASE_Spegel2013(end, 25));
        Spegel_6_pyr = (END_Spegel2013(6,14)+END_Spegel2013(6,25)) / (BASE_Spegel2013(end, 14)+BASE_Spegel2013(end, 25));
        Spegel_10_pyr = (END_Spegel2013(10,14)+END_Spegel2013(10,25)) / (BASE_Spegel2013(end, 14)+BASE_Spegel2013(end, 25));
        Spegel_15_pyr = (END_Spegel2013(end,14)+END_Spegel2013(end,25)) / (BASE_Spegel2013(end, 14)+BASE_Spegel2013(end,25));
        
        Spegel_3_r5p = END_Spegel2013(3,20) / BASE_Spegel2013(end, 20);
        Spegel_6_r5p = END_Spegel2013(6,20) / BASE_Spegel2013(end, 20);
        Spegel_10_r5p = END_Spegel2013(10,20) / BASE_Spegel2013(end, 20);
        Spegel_15_r5p = END_Spegel2013(end,20) / BASE_Spegel2013(end, 20);
        
        Spegel_3_suc = END_Spegel2013(3,31) / BASE_Spegel2013(end, 31);
        Spegel_6_suc = END_Spegel2013(6,31) / BASE_Spegel2013(end, 31);
        Spegel_10_suc = END_Spegel2013(10,31) / BASE_Spegel2013(end, 31);
        Spegel_15_suc = END_Spegel2013(end,31) / BASE_Spegel2013(end, 31);
        
        
        Spegel_3_G6P = END_Spegel2013(3,3) / BASE_Spegel2013(3,3);
        Spegel_6_G6P = END_Spegel2013(6,3) / BASE_Spegel2013(6,3);
        Spegel_10_G6P = END_Spegel2013(10,3) / BASE_Spegel2013(10,3);
        Spegel_15_G6P = END_Spegel2013(15,3) / BASE_Spegel2013(15,3);
        
        Spegel_3_GLC = END_Spegel2013(3,1) / BASE_Spegel2013(3,1);
        Spegel_6_GLC = END_Spegel2013(6,1) / BASE_Spegel2013(6,1);
        Spegel_10_GLC = END_Spegel2013(10,1) / BASE_Spegel2013(10,1);
        Spegel_15_GLC = END_Spegel2013(15,1) / BASE_Spegel2013(15,1);
        
        Spegel_3_glut = (END_Spegel2013(3,36)+END_Spegel2013(3,38)) / (BASE_Spegel2013(end, 36)+BASE_Spegel2013(end, 36));
        Spegel_6_glut = (END_Spegel2013(6,36)+END_Spegel2013(6,38)) / (BASE_Spegel2013(end, 36)+BASE_Spegel2013(end, 36));
        Spegel_10_glut = (END_Spegel2013(10,36)+END_Spegel2013(10,38)) / (BASE_Spegel2013(end, 36)+BASE_Spegel2013(end, 36));
        Spegel_15_glut = (END_Spegel2013(end,36)+END_Spegel2013(end,38)) / (BASE_Spegel2013(end, 36)+BASE_Spegel2013(end,36));
        
    catch
        err = 1;
        fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
        Spegel_3_2pg = 0;
        Spegel_6_2pg = 0;
        Spegel_10_2pg = 0;
        Spegel_15_2pg = 0;
        
        Spegel_3_3pg = 0;
        Spegel_6_3pg = 0;
        Spegel_10_3pg = 0;
        Spegel_15_3pg = 0;
        
        Spegel_3_akg = 0;
        Spegel_6_akg = 0;
        Spegel_10_akg = 0;
        Spegel_15_akg = 0;
        
        Spegel_3_ala =0;
        Spegel_6_ala = 0;
        Spegel_10_ala = 0;
        Spegel_15_ala = 0;
        
        Spegel_3_asp = 0;
        Spegel_6_asp = 0;
        Spegel_10_asp = 0;
        Spegel_15_asp = 0;
        
        Spegel_3_cit = 0;
        Spegel_6_cit = 0;
        Spegel_10_cit = 0;
        Spegel_15_cit = 0;
        
        Spegel_3_fum = 0;
        Spegel_6_fum = 0;
        Spegel_10_fum =0;
        Spegel_15_fum = 0;
        
        Spegel_3_g3p = 0;
        Spegel_6_g3p = 0;
        Spegel_10_g3p = 0;
        Spegel_15_g3p =0;
        
        Spegel_3_lac = 0;
        Spegel_6_lac = 0;
        Spegel_10_lac = 0;
        Spegel_15_lac = 0;
        
        Spegel_3_mal = 0;
        Spegel_6_mal = 0;
        Spegel_10_mal = 0;
        Spegel_15_mal = 0;
        
        Spegel_3_pep = 0;
        Spegel_6_pep = 0;
        Spegel_10_pep = 0;
        Spegel_15_pep = 0;
        
        Spegel_3_pyr = 0;
        Spegel_6_pyr = 0;
        Spegel_10_pyr = 0;
        Spegel_15_pyr = 0;
        
        Spegel_3_r5p = 0;
        Spegel_6_r5p = 0;
        Spegel_10_r5p = 0;
        Spegel_15_r5p = 0;
        
        Spegel_3_suc = 0;
        Spegel_6_suc = 0;
        Spegel_10_suc = 0;
        Spegel_15_suc = 0;
        
        
        Spegel_3_G6P = 0;
        Spegel_6_G6P = 0;
        Spegel_10_G6P = 0;
        Spegel_15_G6P = 0;
        
        Spegel_3_GLC = 0;
        Spegel_6_GLC = 0;
        Spegel_10_GLC = 0;
        Spegel_15_GLC = 0;
        
        Spegel_3_glut = 0;
        Spegel_6_glut = 0;
        Spegel_10_glut =0;
        Spegel_15_glut = 0;
    end
    
    
    %% Graham 30 minute
    try
        tspan = [0:30];
        % Afirst = A;
        initvalue = A;
        params(71,1) = 1e-6;
        [~, BASE] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
        %% Low Glucose no Ex, extGluc = 2.8
        tspan = [0:30];
        %initvalue = (BASE(end,:))';
        %Change this for L,M, H
        params(71,1) = 2.8;
        initvalue = (BASE(end,:))';
        [time, END_GrahamL] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
        
        tspan = [0:30];
        %initvalue = (BASE(end,:))';
        %Change this for L,M, H
        params(71,1) = 16.7;
        initvalue = (BASE(end,:))';
        [time, END_GrahamH] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
        Graham_2PG = END_GrahamH(end,12) / END_GrahamL(end, 12);
        Graham_3PG = END_GrahamH(end,11) / END_GrahamL(end,11);
        Graham_AKG = (END_GrahamH(end, 41) + END_GrahamH(end, 29)) / (END_GrahamL(end, 41) + END_GrahamL(end, 29));
        Graham_ALA = END_GrahamH(end, 50) / END_GrahamL(end, 50);
        Graham_ASP = (END_GrahamH(end,35) + END_GrahamH(end,37)) / (END_GrahamL(end,35) + END_GrahamL(end,37));
        Graham_FUM = END_GrahamH(20,32) / END_GrahamL(20,32);
        Graham_G6P = END_GrahamH(end,3) / END_GrahamL(end,3);
        Graham_GLC = END_GrahamH(end,1) / END_GrahamL(end,1);
        Graham_glut= (END_GrahamH(end, 36) + END_GrahamH(end, 38)) / (END_GrahamL(end, 36) + END_GrahamL(end, 38));
        Graham_G3P = END_GrahamH(end, 8) / END_GrahamL(end, 8);
        Graham_LAC = END_GrahamH(end,15) / END_GrahamL(end,15);
        Graham_MAL = (END_GrahamH(end,33) + END_GrahamH(end,40)) / (END_GrahamL(end,33) + END_GrahamL(end,40));
        Graham_SUC = END_GrahamH(end,31) / END_GrahamL(end,31);
        Graham_CIT = (END_GrahamH(end,27)+END_GrahamH(end,42)) / (END_GrahamL(end, 27)+END_GrahamL(end, 42));
        Graham_PEP = END_GrahamH(end,13) / END_GrahamL(end, 13);
        Graham_PYR = (END_GrahamH(end,14)+END_GrahamH(end,25)) / (END_GrahamL(end, 14)+END_GrahamL(end,25));
        Graham_R5P = END_GrahamH(end,20) / END_GrahamH(end, 20);
    catch
        err = 1;
        fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
        Graham_2PG = 0;
        Graham_3PG = 0;
        Graham_AKG = 0;
        Graham_ALA = 0;
        Graham_ASP = 0;
        Graham_FUM = 0;
        Graham_G6P = 0;
        Graham_GLC = 0;
        Graham_glut= 0;
        Graham_G3P = 0;
        Graham_LAC =0;
        Graham_MAL = 0;
        Graham_SUC = 0;
        Graham_CIT = 0;
        Graham_PEP = 0;
        Graham_PYR = 0;
        Graham_R5P = 0;
    end
    %% Malmgren-60mins
    
    try
        tspan = [0:120];
        initvalue = A;
        params(71,1) = 2.8;
        [~, BASE] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
        
        tspan = [0:60];
        params(71,1) = 2.8;
        initvalue = (BASE(end,:))';
        [time, END_MalmgrenL] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
        
        params(71,1) = 16.7;
        initvalue = (BASE(end,:))';
        [time, END_MalmgrenH] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
        
        % Metformin 60min
        params(71,1) = 16.7;
        initvalue = (BASE(end,:))';
        [time, END_MalmgrenH_Met] = ode15s(@INS1Epathway_Metformin,tspan,initvalue,options,params);
        
        Malmgren_AKG = END_MalmgrenH(end, 41) + END_MalmgrenH(end, 29) / END_MalmgrenL(end, 41) + END_MalmgrenL(end, 29);
        Malmgren_ALA = END_MalmgrenH(end, 50) / END_MalmgrenL(end, 50);
        Malmgren_ASP = END_MalmgrenH(end,35) + END_MalmgrenH(end,37) / END_MalmgrenL(end,35) + END_MalmgrenL(end,37);
        Malmgren_CIT = (END_MalmgrenH(end,27)+END_MalmgrenH(end,42)) / (END_MalmgrenL(end, 27)+END_MalmgrenL(end, 42));
        Malmgren_FUM = END_MalmgrenH(end,32) / END_MalmgrenL(end,32);
        Malmgren_G6P = END_MalmgrenH(end,3) / END_MalmgrenL(end,3);
        Malmgren_GLC = END_MalmgrenH(end,1) / END_MalmgrenL(end,1);
        Malmgren_glut= (END_MalmgrenH(end, 36) + END_MalmgrenH(end, 38)) / (END_MalmgrenL(end, 36) + END_MalmgrenL(end, 38));
        Malmgren_G3P = END_MalmgrenH(end, 8) / END_MalmgrenL(end, 8);
        Malmgren_ICIT= (END_MalmgrenH(end,28)+END_MalmgrenH(end,51)) / (END_MalmgrenL(end, 28)+END_MalmgrenL(end, 51));
        Malmgren_LAC = END_MalmgrenH(end,15) / END_MalmgrenL(end,15);
        Malmgren_MAL = (END_MalmgrenH(end,33) + END_MalmgrenH(end,40)) / (END_MalmgrenL(end,33) + END_MalmgrenL(end,40));
        Malmgren_PYR = (END_MalmgrenH(end,14)+END_MalmgrenH(end,25)) / (END_MalmgrenL(end, 14)+END_MalmgrenL(end, 25));
        Malmgren_SUC = END_MalmgrenH(end,31) / END_MalmgrenL(end,31);
        Malmgren_2PG = END_MalmgrenH(end,12) / END_MalmgrenL(end,12);
        Malmgren_3PG = END_MalmgrenH(end,11) / END_MalmgrenL(end,11);
        Malmgren_PEP = END_MalmgrenH(end,13) / END_MalmgrenL(end,13);
        Malmgren_R5P = END_MalmgrenH(end,20) / END_MalmgrenL(end,20);
        %% Goehring 45min
        Goehring_2PG = END_MalmgrenH(45,12) / END_MalmgrenL(45,12);
        Goehring_3PG = END_MalmgrenH(45,11) / END_MalmgrenL(45,11);
        Goehring_3PG_60 = END_MalmgrenH(end,11) / END_MalmgrenL(end,11);
        Goehring_AKG = END_MalmgrenH(45, 41) + END_MalmgrenH(45, 29) / END_MalmgrenL(45, 41) + END_MalmgrenL(45, 29);
        Goehring_CIT = (END_MalmgrenH(45,27)+END_MalmgrenH(45,42)) / (END_MalmgrenL(45, 27)+END_MalmgrenL(45, 42));
        Goehring_FUM = END_MalmgrenH(45,32) / END_MalmgrenL(45,32);
        Goehring_G6P = END_MalmgrenH(45,3) / END_MalmgrenL(45,3);
        Goehring_GLC = END_MalmgrenH(45,1) / END_MalmgrenL(45,1);
        Goehring_LAC = END_MalmgrenH(45,15) / END_MalmgrenL(45,15);
        Goehring_MAL = (END_MalmgrenH(45,33) + END_MalmgrenH(45,40)) / (END_MalmgrenL(45,33) + END_MalmgrenL(45,40));
        Goehring_PEP = END_MalmgrenH(45,13) / END_MalmgrenL(45,13);
        Goehring_PYR = (END_MalmgrenH(45,14)+END_MalmgrenH(45,25)) / (END_MalmgrenL(45, 14)+END_MalmgrenL(45, 25));
        Goehring_R5P = END_MalmgrenH(45,20) / END_MalmgrenL(45,20);
        Goehring_SUC = END_MalmgrenH(45,31) / END_MalmgrenL(45,31);
        Goehring_ALA = END_MalmgrenH(45, 50) / END_MalmgrenL(45, 50);
        Goehring_glut= (END_MalmgrenH(45, 36) + END_MalmgrenH(45, 38)) / (END_MalmgrenL(45, 36) + END_MalmgrenL(45, 38));
        Goehring_ASP = END_MalmgrenH(45,35) + END_MalmgrenH(45,37) / END_MalmgrenL(45,35) + END_MalmgrenL(45,37);
        
    catch
        err = 1;
        fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
        Malmgren_AKG = 0;
        Malmgren_ALA = 0;
        Malmgren_ASP = 0;
        Malmgren_CIT = 0;
        Malmgren_FUM = 0;
        Malmgren_G6P = 0;
        Malmgren_GLC = 0;
        Malmgren_glut= 0;
        Malmgren_G3P = 0;
        Malmgren_LAC = 0;
        Malmgren_MAL = 0;
        Malmgren_PYR = 0;
        Malmgren_SUC = 0;
        Malmgren_2PG = 0;
        Malmgren_3PG = 0;
        Malmgren_PEP = 0;
        Malmgren_R5P = 0;
        %% Goehring 45min
        Goehring_2PG = 0;
        Goehring_3PG = 0;
        Goehring_3PG_60 = 0;
        Goehring_AKG = 0;
        Goehring_CIT = 0;
        Goehring_FUM = 0;
        Goehring_G6P =0;
        Goehring_GLC = 0;
        Goehring_LAC = 0;
        Goehring_MAL = 0;
        Goehring_PEP = 0;
        Goehring_PYR = 0;
        Goehring_R5P = 0;
        Goehring_SUC = 0;
        Goehring_ALA =0;
        Goehring_glut= 0;
        Goehring_ASP = 0;
    end
%     observables(:,i) = [Spegel_3_2pg
%         Spegel_6_2pg
%         Spegel_10_2pg
%         Spegel_15_2pg
%         Graham_2PG
%         Goehring_2PG
%         Malmgren_2PG
%         Spegel_3_3pg
%         Spegel_6_3pg
%         Spegel_10_3pg
%         Spegel_15_3pg
%         Graham_3PG
%         Goehring_3PG
%         Malmgren_3PG
%         Spegel_3_akg
%         Spegel_6_akg
%         Spegel_10_akg
%         Spegel_15_akg
%         Graham_AKG
%         Goehring_AKG
%         Malmgren_AKG
%         Spegel_3_ala
%         Spegel_6_ala
%         Spegel_10_ala
%         Spegel_15_ala
%         Graham_ALA
%         Goehring_ALA
%         Malmgren_ALA
%         Spegel_3_asp
%         Spegel_6_asp
%         Spegel_10_asp
%         Spegel_15_asp
%         Graham_ASP
%         Goehring_ASP
%         Malmgren_ASP
%         Spegel_3_cit
%         Spegel_6_cit
%         Spegel_10_cit
%         Spegel_15_cit
%         Graham_CIT
%         Goehring_CIT
%         Malmgren_CIT
%         Spegel_3_fum
%         Spegel_6_fum
%         Spegel_10_fum
%         Spegel_15_fum
%         Graham_FUM
%         Goehring_FUM
%         Malmgren_FUM
%         Spegel_3_g3p
%         Spegel_6_g3p
%         Spegel_10_g3p
%         Spegel_15_g3p
%         Graham_G3P
%         Goehring_3PG
%         Malmgren_G3P
%         Spegel_3_lac
%         Spegel_6_lac
%         Spegel_10_lac
%         Spegel_15_lac
%         Graham_LAC
%         Goehring_LAC
%         Malmgren_LAC
%         Spegel_3_mal
%         Spegel_6_mal
%         Spegel_10_mal
%         Spegel_15_mal
%         Graham_MAL
%         Goehring_MAL
%         Malmgren_MAL
%         Spegel_3_pep
%         Spegel_6_pep
%         Spegel_10_pep
%         Spegel_15_pep
%         Graham_PEP
%         Goehring_PEP
%         Malmgren_PEP
%         Spegel_3_pyr
%         Spegel_6_pyr
%         Spegel_10_pyr
%         Spegel_15_pyr
%         Graham_PYR
%         Goehring_PYR
%         Malmgren_PYR
%         Spegel_3_r5p
%         Spegel_6_r5p
%         Spegel_10_r5p
%         Spegel_15_r5p
%         Graham_R5P
%         Goehring_R5P
%         Malmgren_R5P
%         Spegel_3_suc
%         Spegel_6_suc
%         Spegel_10_suc
%         Spegel_15_suc
%         Graham_SUC
%         Goehring_SUC
%         Malmgren_SUC
%         Spegel_3_G6P
%         Spegel_6_G6P
%         Spegel_10_G6P
%         Spegel_15_G6P
%         Graham_G6P
%         Goehring_G6P
%         Malmgren_G6P
%         Spegel_3_GLC
%         Spegel_6_GLC
%         Spegel_10_GLC
%         Spegel_15_GLC
%         Graham_GLC
%         Goehring_GLC
%         Malmgren_GLC
%         Spegel_3_glut
%         Spegel_6_glut
%         Spegel_10_glut
%         Spegel_15_glut
%         Graham_glut
%         Goehring_glut
%         Malmgren_glut]';
        Spegel_3 = (END_Spegel2013(3,:) ./ BASE_Spegel2013(end, :))';
        Spegel_6 = (END_Spegel2013(6,:) ./ BASE_Spegel2013(end, :))';
        Spegel_10 = (END_Spegel2013(10,:) ./ BASE_Spegel2013(end, :))';
        Spegel_15 = (END_Spegel2013(end,:) ./ BASE_Spegel2013(end, :))';
        Malmgren = (END_MalmgrenH(end,:) ./ END_MalmgrenL(end,:))';
        MalmgrenMet = (END_MalmgrenH_Met(end,:) ./ END_MalmgrenH(end,:))';
    
    observables = horzcat(Spegel_3, Spegel_6, Spegel_10,Spegel_15,Malmgren, MalmgrenMet);  
end