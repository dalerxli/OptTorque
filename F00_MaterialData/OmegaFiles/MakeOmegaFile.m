%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Matlab code to Create an OmegaFile                      %%%%%%%
%%%%%%% Yoonkyung Eunnie Lee 2015-03-11                         %%%%%%%
%%% function MakeOmegaFile(lambdamin, lambdamax, numpts, OmegaFileName)
%%% lambdamin, lambdamax: wavelength in [nm]
%%% numpts: number of data points 
%%% OmegaFileName: ex. 'OmegaFile_300' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeOmegaFile(lambdamin, lambdamax, numpts, OmegaFileName)

% Omegafile : omega is in the order of 2*pi*10^(-6)*/lambda.
lambda = linspace(lambdamin, lambdamax, numpts) ;
omega = 2*pi*1000./lambda;
omega = omega' ;

dlmwrite(OmegaFileName, omega,'delimiter','\t','precision','%6.5f'); 
end