%% Setup
function [aveSpegel3,aveSpegel6, aveSpegel10, aveSpegel15, aveGoehring45, aveMalmgren60, allAve, Spegel3, Spegel6, Spegel10, Spegel15, Goehring45, Malmgren60 ] = GetAveRxnRate(params)
A = [0.542125648135445;3.76168895009916;0.454448540472248;0.706990028881560;0.0156565460978362;0.478582126811771;4.39031164654454;0.168038706461350;4.21968660787895;0.0560506065370617;0.397377606983681;0.583587181302951;0.716708019845057;0.000154147902632731;0.136501997249423;0.132939511019276;0.00456400046649098;0.000733943694403398;0.00134582377890930;0.00483413989498314;0.0379963772949342;0.00192408120909374;4.23518730127469;2.45107434802339e-06;0.000771217660846397;2.09292143332501;0.288669565535384;0.00779378981985000;0.224895886981259;12.6827767477872;6.03285037234089e-05;6.77882129684873e-05;0.000309121127035145;1.74318833008063e-06;1.45680249381982e-05;2.94658145921870;3.12055381157148e-05;0.797796208874077;2.65472943228149e-07;2.85338922615892e-05;0.0240841285605586;1.05706088817784e-05;3.31601148362778;0.0170013921210288;0.00150069872532053;2.47590966892564;0.0300000000000000;2.15537701634457e-07;4.94999978446232;3.43977996188819;1.51575084594245e-06;0.0400301566013678;1.72526016292169;0.0300000000000000;1.03960000000000;0.0480000000000000];
options = odeset('RelTol',1e-12,...
    'AbsTol',1e-19,...
    'NormControl','on',...
    'Stats','off',...
    'BDF','off',...
    'NonNegative',1:56,...
    'MaxOrder', 5);

%% Spegel 2013
tspan = [0:136];
try
    % Afirst = A;
    initvalue = A;
    params(71,1) = 2.8;
    [~, BASE_Spegel2013] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
    
    tspan = [0:15];
    initvalue = (BASE_Spegel2013(120,:))';
    params(71,1) = 16.7;
    [time, END_Spegel2013] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
    Spegel3 = [];
    Spegel6 = [];
    Spegel10 = [];
    Spegel15 = [];
    
    for i = 1:3
        y = BASE_Spegel2013(i+120,:);
        [reactions_LG] = getRXN_Velocity_INS1E(y,params);
        y = END_Spegel2013(i,:);
        [reactions_HG] = getRXN_Velocity_INS1E(y,params);
        Spegel3(:,i) = reactions_HG ./ reactions_LG;
    end
    Spegel3(isnan(Spegel3))=0.0001;
    aveSpegel3 = mean(Spegel3, 2, 'omitnan');
    for i = 1:6
        y = BASE_Spegel2013(i+120,:);
        [reactions_LG] = getRXN_Velocity_INS1E(y,params);
        y = END_Spegel2013(i,:);
        [reactions_HG] = getRXN_Velocity_INS1E(y,params);
        Spegel6(:,i) = reactions_HG ./ reactions_LG;
    end
    aveSpegel6 = mean(Spegel6, 2, 'omitnan');
    Spegel6(isnan(Spegel6))=0.0001;
    for i = 1:10
        y = BASE_Spegel2013(i+120,:);
        [reactions_LG] = getRXN_Velocity_INS1E(y,params);
        y = END_Spegel2013(i,:);
        [reactions_HG] = getRXN_Velocity_INS1E(y,params);
        Spegel10(:,i) = reactions_HG ./ reactions_LG;
    end
    aveSpegel10 = mean(Spegel10, 2, 'omitnan');
    Spegel10(isnan(Spegel10))=0.0001;
    
    for i = 1:15
        y = BASE_Spegel2013(i+120,:);
        [reactions_LG] = getRXN_Velocity_INS1E(y,params);
        y = END_Spegel2013(i,:);
        [reactions_HG] = getRXN_Velocity_INS1E(y,params);
        Spegel15(:,i) = reactions_HG ./ reactions_LG;
    end
    aveSpegel15 = mean(Spegel15, 2, 'omitnan');
    Spegel15(isnan(Spegel15))=0.0001;
    
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    Spegel3 =zeros(61,3);
    Spegel6 = zeros(61,6);
    Spegel10 = zeros(61,10);
    Spegel15 = zeros(61,15);
    aveSpegel3 = zeros(61,1);
    aveSpegel6 = zeros(61,1);
    aveSpegel10 =zeros(61,1);
    aveSpegel15 = zeros(61,1);
end

%% Malmgren-60mins
try
    tspan = [0:120];
    initvalue = A;
    params(71,1) = 2.8;
    [~, BASE] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
    tspan = [0:65];
    params(71,1) = 2.8;
    initvalue = (BASE(end,:))';
    [time, END_MalmgrenL] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
    params(71,1) = 16.7;
    initvalue = (BASE(end,:))';
    [time, END_MalmgrenH] = ode15s(@INS1Epathway,tspan,initvalue,options,params);
    params(71,1) = 16.7;
    initvalue = (BASE(end,:))';
    [time, END_MalmgrenH_Met] = ode15s(@INS1Epathway_Metformin,tspan,initvalue,options,params);
    Goehring45 = [];
    Malmgren60 = [];
    for i = 1:45
        y = END_MalmgrenL(i,:);
        [reactions_LG] = getRXN_Velocity_INS1E(y,params);
        y = END_MalmgrenH(i,:);
        [reactions_HG] = getRXN_Velocity_INS1E_Metformin(y,params);
        Goehring45(:,i) = reactions_HG ./ reactions_LG;
    end
    
    
    aveGoehring45 = mean(Goehring45, 2, 'omitnan');
    Goehring45(isnan(Goehring45))=0.0001;
    
    for i = 1:60
        y = END_MalmgrenH(i,:);
        [reactions_LG] = getRXN_Velocity_INS1E(y,params);
        y = END_MalmgrenH_Met(i,:);
        [reactions_HG] = getRXN_Velocity_INS1E_Metformin(y,params);
        Malmgren60(:,i) = reactions_HG ./ reactions_LG;
    end
    aveMalmgren60 = mean(Malmgren60, 2, 'omitnan');
    Malmgren60(isnan(Malmgren60))=0.0001;
    
    
    allAve = horzcat(aveSpegel3, aveSpegel6, aveSpegel10, aveSpegel15, ...
        aveMalmgren60);
    allAve(isnan(allAve))=0.0001;
    
    
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    Goehring45 = zeros(61,45);
    Malmgren60 = zeros(61,60);
    aveGoehring45 =zeros(61,1);
    aveMalmgren60 = zeros(61,1);
    allAve = horzcat(aveSpegel3, aveSpegel6, aveSpegel10, aveSpegel15, ...
        aveMalmgren60);
    
end
end
