%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The 'shocks' block specifies the non zero elements of the covariance 
% matrix of the shocks of exogenous variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2013 Benjamin K. Johannsen, Lawrence J. Christiano
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see http://www.gnu.org/licenses/.

shocks; 
var e_xp;       stderr stdxp_p;
var e_lambdaf;  stderr stdlambdaf_p;
var e_pitarget; stderr stdpitarget_p;
var e_muup;     stderr stdmuup_p;
var e_g;        stderr stdg_p;
var e_s;  stderr stdmuzstar_p;
@# if cee == 0
var e_gamma;    stderr stdgamma_p;
@# endif
@# if cee == 0
var e_sigma;    stderr 1-@{stopunant};
var e_sigma2;    stderr 1-@{stopunant};
var e_sigma3;    stderr 1-@{stopunant};
@# endif
var e_zetac;    stderr stdzetac_p;
@# if cee == 0
@#if Spread1_in_financial_data
var e_term;     stderr stdterm_p;
@#else
var e_term;     stderr 0;
@#endif
var e_xi8;      stderr 1-@{stopsignal};
var e_xi7;      stderr 1-@{stopsignal};
var e_xi6;      stderr 1-@{stopsignal};
var e_xi5;      stderr 1-@{stopsignal};
var e_xi4;      stderr 1-@{stopsignal};
var e_xi3;      stderr 1-@{stopsignal};
var e_xi2;      stderr 1-@{stopsignal};
var e_xi1;      stderr 1-@{stopsignal};
var e_xi82;      stderr 1-@{stopsignal};
var e_xi72;      stderr 1-@{stopsignal};
var e_xi62;      stderr 1-@{stopsignal};
var e_xi52;      stderr 1-@{stopsignal};
var e_xi42;      stderr 1-@{stopsignal};
var e_xi32;      stderr 1-@{stopsignal};
var e_xi22;      stderr 1-@{stopsignal};
var e_xi12;      stderr 1-@{stopsignal};
var e_xi83;      stderr 1-@{stopsignal};
var e_xi73;      stderr 1-@{stopsignal};
var e_xi63;      stderr 1-@{stopsignal};
var e_xi53;      stderr 1-@{stopsignal};
var e_xi43;      stderr 1-@{stopsignal};
var e_xi33;      stderr 1-@{stopsignal};
var e_xi23;      stderr 1-@{stopsignal};
var e_xi13;      stderr 1-@{stopsignal};
@# endif
end;