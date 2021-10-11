%% mnb_Simple_PSO_ObjFunc
% The objective function used by mnb_Simple_PSO_driver(). This function is
% called from within mnb_PSOFit().
%
% Usage: fx = mnb_Simple_PSO_ObjFunc(pos,ObjFuncData)
%

%% Required Input Arguments (in order)
%
% * pos : numeric vector : Positions of the particle (the values of the fit
%       parameters).
%
% * ObjFuncData : any variable type : The argument is simply passed-through
%       from mnb_PSOFit() to this objective function. It can be used, for
%       example, to store the experimental data, a SimBiology model, etc.
%
%% Optional Input Arguments (passed via Param-Value pairs)
%
%
%% Function Output Arguments
%
% * fx : numeric : calculated fitness value
%
%% Dependencies
%
%
%% Limitations
%
%
%% Version Information
% * Last changed date: $Date: 2015-03-04 10:49:38 -0500 (Wed, 04 Mar 2015) $
% * Last changed by: $Author: jkearns $
% * Last changed revision: $Revision: 1560 $
%
% Copyright 2012-2015 Merrimack Pharmaceuticals, Inc.

function [fx] = ODE_ObjFunc_CollectedTimeCourse(pos, ObjFuncData)
%[fx] = ODE_ObjFunc_MG(pos,ObjFuncData)
% This example uses two matrices in the ObjFuncData input argument
%
%  ObjFuncData.expdata - 2 column matrix in which the first column is
%                   the drug concentration and the 2nd is the mean
%                   value of the response.
%
%  ObjFuncData.Stdev - 2 column matrix in which the first column is
%                   the drug concentration and the 2nd is the standard
%                   deviation of the response.


% Compute the predicted response value for each drug concentration
% compdata = (100*(1-ObjFuncData.expdata(:,1).^pos(1)./...
% (pos(2)+ObjFuncData.expdata(:,1).^pos(1))))*...
% (100-pos(3))/100+pos(3);
try
    [~, observables_out] = simulation_protocol_fitting_CollectedTimeCourse(pos);
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    %observables_out = zeros(1,24);
    observables_out = zeros(81,1);
end

compdata = (observables_out);
% Evaluate the fitness
% fx = sum(((ObjFuncData(:,1) - compdata(:,1))./ObjFuncData(:,1)).^2);

fx = 0;
for i = 1:length(ObjFuncData) 
fx = fx + ((1/ObjFuncData(i,2))*((ObjFuncData(i,1) - compdata(i,1))^2));
end
% normalizing data by magnitude and stdev
%fx = sum((((ObjFuncData(:,1) - compdata(:,1))./ObjFuncData(:,1)).^2) .* (1./(ObjFuncData(:,2))));
% Normalizing by variance 
%fx = sum(((ObjFuncData(:,1) - compdata(:,1)).^2) .* (1./(ObjFuncData(:,2).^2)));

end
