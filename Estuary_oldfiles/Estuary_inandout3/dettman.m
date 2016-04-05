function [N1]=dettman(tc1,R,N0,B)
% dettman- this is a placeholder for the process that defines the stock of N
% in the estuary in a year (N1)

% tc1 is the total loading in a year
% R the % carryover from last year
% N0 N stock last year
% B is the background loading not from the controllable sources


N1 = (R*N0 + tc1)+B;
 
end
