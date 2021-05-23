function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = cmr.static_resid_tt(T, y, x, params);
end
residual = zeros(349, 1);
lhs = y(28);
rhs = T(41);
residual(1) = lhs - rhs;
lhs = y(5);
rhs = T(18)*y(55)*y(17)+y(5)*params(71)*T(5)*params(8);
residual(2) = lhs - rhs;
lhs = y(5)*T(7);
rhs = y(37)*T(18)*y(17)*y(16)*y(55)+y(5)*T(7)*T(43);
residual(3) = lhs - rhs;
lhs = y(6);
rhs = y(10)*T(16)*y(55)*y(17)*(1-params(66))/params(22)+y(6)*T(51);
residual(4) = lhs - rhs;
lhs = T(30);
rhs = y(55)*T(52)*params(73)+T(30)*T(53);
residual(5) = lhs - rhs;
lhs = y(54);
rhs = T(54)^(1/T(49));
residual(6) = lhs - rhs;
lhs = y(34);
rhs = exp(params(49)*(y(51)-1))*(y(34))*params(67);
residual(7) = lhs - rhs;
lhs = y(34);
rhs = y(37)*y(4)*params(4)*T(55);
residual(8) = lhs - rhs;
lhs = y(37);
rhs = T(58)/y(4);
residual(9) = lhs - rhs;
lhs = T(18);
rhs = params(9)*(1-y(8))*(y(20)-params(70))/y(8)+y(15)*y(29)*(1+y(35))*(T(12)-T(13))/(y(24)*y(19))+y(1)+y(7)+y(12)/y(18)+y(15)*T(59)/(y(19)*params(69));
residual(10) = lhs - rhs;
lhs = y(15);
rhs = y(15)*(1-params(11))/(y(19)*params(69))+y(12)*(1-(exp(T(31)*(params(69)*y(19)*y(56)-params(69)*params(25)))+exp((params(69)*y(19)*y(56)-params(69)*params(25))*(-T(31)))-2));
residual(11) = lhs - rhs;
lhs = 0;
rhs = y(17)*y(55)*params(8)/(y(24)*y(19))*(1+(1-params(64))*y(30))-y(55)*y(17);
residual(12) = lhs - rhs;
lhs = y(17)*y(55)*(1+params(63));
rhs = y(19)*y(55)/(y(19)*y(1)-y(1)*params(7))-y(55)*params(8)*params(7)/(y(19)*y(1)-y(1)*params(7));
residual(13) = lhs - rhs;
lhs = 0;
rhs = (1+y(35))*(1-T(12))/(1+y(30))+T(61)*T(62);
residual(14) = lhs - rhs;
lhs = 1+y(35);
rhs = params(11)*params(65)+y(24)*(y(29)*(1-params(11))+(1-params(65))*(y(34)*y(51)-T(59)))/(y(29)*params(69));
residual(15) = lhs - rhs;
lhs = 0;
rhs = y(17)*(-y(55))/y(18)+y(29)*y(55)*y(17)*T(64)+T(65)*T(66)/(y(19)*params(69));
residual(16) = lhs - rhs;
lhs = log(y(30)/params(33));
rhs = log(y(30)/params(33))*params(45)+1/params(33)*(1-params(45))*params(27)*log(y(26)/params(27))+T(67)*(log(y(24))-log(y(26)))+params(25)*(1-params(45))*1/(params(33)*4)*params(3)*log(T(69))+log(T(69))*params(25)*1/params(33)*(1-params(45))*params(1)-T(70)*(params(10)*log(y(1)/params(10))+params(17)*log(y(12)/params(17))-params(17)*log(y(18))+params(15)*log(y(7)/params(15)))/T(68)+1/(params(33)*400)*x(18);
residual(17) = lhs - rhs;
residual(18) = 1+T(13)*(1+y(35))*y(15)*y(29)/(y(20)*(1+y(30)))-y(15)*y(29)/y(20);
lhs = y(20);
rhs = y(20)*y(8)*(1+y(30))/(y(24)*y(19))+params(70)+y(29)*y(15)*T(71)*(y(35)-y(30)-(1+y(35))*(T(12)-T(13)));
residual(19) = lhs - rhs;
lhs = y(50);
rhs = T(72)*T(74);
residual(20) = lhs - rhs;
lhs = y(55)*y(17);
rhs = T(119);
residual(21) = lhs - rhs;
lhs = y(55)*y(17);
rhs = T(141);
residual(22) = lhs - rhs;
lhs = y(2);
rhs = T(69);
residual(23) = lhs - rhs;
lhs = y(3);
rhs = T(69);
residual(24) = lhs - rhs;
lhs = y(9);
rhs = T(69);
residual(25) = lhs - rhs;
lhs = y(11);
rhs = y(10)/(y(10));
residual(26) = lhs - rhs;
lhs = y(13);
rhs = y(24)/params(27);
residual(27) = lhs - rhs;
lhs = y(14);
rhs = T(69);
residual(28) = lhs - rhs;
lhs = y(21);
rhs = T(69);
residual(29) = lhs - rhs;
lhs = y(27);
rhs = T(76);
residual(30) = lhs - rhs;
lhs = y(25);
rhs = 1;
residual(31) = lhs - rhs;
lhs = y(31);
rhs = exp(y(30)-params(33));
residual(32) = lhs - rhs;
lhs = y(32);
rhs = (1+y(30))/y(24)/((1+params(33))/params(27));
residual(33) = lhs - rhs;
lhs = y(48);
rhs = 1+y(36)-y(30);
residual(34) = lhs - rhs;
lhs = y(52);
rhs = T(69);
residual(35) = lhs - rhs;
lhs = log(y(4)/params(12));
rhs = log(y(4)/params(12))*params(34)+x(1);
residual(36) = lhs - rhs;
lhs = T(77);
rhs = x(2)+T(77)*params(35);
residual(37) = lhs - rhs;
lhs = log(y(8)/params(16));
rhs = log(y(8)/params(16))*params(36)+x(3);
residual(38) = lhs - rhs;
lhs = log(y(16)/params(21));
rhs = log(y(16)/params(21))*params(37)+x(4);
residual(39) = lhs - rhs;
lhs = log(y(18)/params(24));
rhs = log(y(18)/params(24))*params(38)+x(5);
residual(40) = lhs - rhs;
lhs = log(T(69));
rhs = log(T(69))*params(39)+x(6);
residual(41) = lhs - rhs;
lhs = log(y(26)/params(29));
rhs = log(y(26)/params(29))*params(40)+x(7);
residual(42) = lhs - rhs;
lhs = log(y(49)/params(68));
rhs = log(y(49)/params(68))*params(42)+x(17);
residual(43) = lhs - rhs;
lhs = log(y(55)/params(74));
rhs = log(y(55)/params(74))*params(43)+x(19);
residual(44) = lhs - rhs;
lhs = log(y(56)/params(75));
rhs = log(y(56)/params(75))*params(44)+x(20);
residual(45) = lhs - rhs;
lhs = log(y(38)/(y(38)));
rhs = log(y(47))+log(y(46))+log(y(45))+log(y(44))+log(y(43))+log(y(42))+log(y(40))+log(y(39))+log(y(38)/(y(38)))*params(41)+log(y(41));
residual(46) = lhs - rhs;
lhs = log(y(47));
rhs = params(62)*x(16);
residual(47) = lhs - rhs;
lhs = log(y(46));
rhs = x(16)*params(62)*params(26)+T(80)*x(15);
residual(48) = lhs - rhs;
lhs = log(y(45));
rhs = x(16)*params(62)*T(78)+x(15)*params(62)*params(26)*T(79)+T(80)*x(14);
residual(49) = lhs - rhs;
lhs = log(y(44));
rhs = x(16)*params(62)*T(81)+x(15)*params(62)*T(78)*T(79)+params(62)*params(26)*T(79)*x(14)+T(80)*x(13);
residual(50) = lhs - rhs;
lhs = log(y(43));
rhs = x(16)*params(62)*T(82)+x(15)*params(62)*T(79)*T(81)+x(14)*params(62)*T(78)*T(79)+params(62)*params(26)*T(79)*x(13)+T(80)*x(12);
residual(51) = lhs - rhs;
lhs = log(y(42));
rhs = x(16)*params(62)*T(83)+x(15)*params(62)*T(79)*T(82)+x(14)*params(62)*T(79)*T(81)+params(62)*T(78)*T(79)*x(13)+params(62)*params(26)*T(79)*x(12)+T(80)*x(11);
residual(52) = lhs - rhs;
lhs = log(y(41));
rhs = x(16)*params(62)*T(84)+x(15)*params(62)*T(79)*T(83)+x(14)*params(62)*T(79)*T(82)+x(13)*params(62)*T(79)*T(81)+params(62)*T(78)*T(79)*x(12)+params(62)*params(26)*T(79)*x(11)+T(80)*x(10);
residual(53) = lhs - rhs;
lhs = log(y(40));
rhs = x(16)*params(62)*T(85)+x(15)*params(62)*T(79)*T(84)+x(14)*params(62)*T(79)*T(83)+x(13)*params(62)*T(79)*T(82)+params(62)*T(79)*T(81)*x(12)+params(62)*T(78)*T(79)*x(11)+params(62)*params(26)*T(79)*x(10)+T(80)*x(9);
residual(54) = lhs - rhs;
lhs = log(y(39));
rhs = x(16)*params(26)^8*params(57)+x(15)*params(57)*T(79)*T(85)+x(14)*T(79)*T(84)*params(57)+x(13)*T(79)*T(83)*params(57)+x(12)*T(79)*T(82)*params(57)+x(11)*T(79)*T(81)*params(57)+x(10)*T(78)*T(79)*params(57)+x(9)*params(26)*T(79)*params(57)+T(79)*params(57)*x(8);
residual(55) = lhs - rhs;
lhs = y(23);
rhs = (y(23));
residual(56) = lhs - rhs;
lhs = y(57);
rhs = T(119);
residual(57) = lhs - rhs;
lhs = y(58);
rhs = T(119);
residual(58) = lhs - rhs;
lhs = y(59);
rhs = T(119);
residual(59) = lhs - rhs;
lhs = y(60);
rhs = T(119);
residual(60) = lhs - rhs;
lhs = y(61);
rhs = T(119);
residual(61) = lhs - rhs;
lhs = y(62);
rhs = T(119);
residual(62) = lhs - rhs;
lhs = y(63);
rhs = T(119);
residual(63) = lhs - rhs;
lhs = y(64);
rhs = T(119);
residual(64) = lhs - rhs;
lhs = y(65);
rhs = T(119);
residual(65) = lhs - rhs;
lhs = y(66);
rhs = T(119);
residual(66) = lhs - rhs;
lhs = y(67);
rhs = T(119);
residual(67) = lhs - rhs;
lhs = y(68);
rhs = T(119);
residual(68) = lhs - rhs;
lhs = y(69);
rhs = T(119);
residual(69) = lhs - rhs;
lhs = y(70);
rhs = T(119);
residual(70) = lhs - rhs;
lhs = y(71);
rhs = T(119);
residual(71) = lhs - rhs;
lhs = y(72);
rhs = T(119);
residual(72) = lhs - rhs;
lhs = y(73);
rhs = T(119);
residual(73) = lhs - rhs;
lhs = y(74);
rhs = T(119);
residual(74) = lhs - rhs;
lhs = y(75);
rhs = T(119);
residual(75) = lhs - rhs;
lhs = y(76);
rhs = T(119);
residual(76) = lhs - rhs;
lhs = y(77);
rhs = T(119);
residual(77) = lhs - rhs;
lhs = y(78);
rhs = T(119);
residual(78) = lhs - rhs;
lhs = y(79);
rhs = T(119);
residual(79) = lhs - rhs;
lhs = y(80);
rhs = T(119);
residual(80) = lhs - rhs;
lhs = y(81);
rhs = T(119);
residual(81) = lhs - rhs;
lhs = y(82);
rhs = T(119);
residual(82) = lhs - rhs;
lhs = y(83);
rhs = T(119);
residual(83) = lhs - rhs;
lhs = y(84);
rhs = T(119);
residual(84) = lhs - rhs;
lhs = y(85);
rhs = T(119);
residual(85) = lhs - rhs;
lhs = y(86);
rhs = T(119);
residual(86) = lhs - rhs;
lhs = y(87);
rhs = T(119);
residual(87) = lhs - rhs;
lhs = y(88);
rhs = T(119);
residual(88) = lhs - rhs;
lhs = y(89);
rhs = T(119);
residual(89) = lhs - rhs;
lhs = y(90);
rhs = T(119);
residual(90) = lhs - rhs;
lhs = y(91);
rhs = T(119);
residual(91) = lhs - rhs;
lhs = y(92);
rhs = T(119);
residual(92) = lhs - rhs;
lhs = y(93);
rhs = T(119);
residual(93) = lhs - rhs;
lhs = y(94);
rhs = T(119);
residual(94) = lhs - rhs;
lhs = y(95);
rhs = T(119);
residual(95) = lhs - rhs;
lhs = y(96);
rhs = T(141);
residual(96) = lhs - rhs;
lhs = y(97);
rhs = T(141);
residual(97) = lhs - rhs;
lhs = y(98);
rhs = T(141);
residual(98) = lhs - rhs;
lhs = y(99);
rhs = T(141);
residual(99) = lhs - rhs;
lhs = y(100);
rhs = T(141);
residual(100) = lhs - rhs;
lhs = y(101);
rhs = T(141);
residual(101) = lhs - rhs;
lhs = y(102);
rhs = T(141);
residual(102) = lhs - rhs;
lhs = y(103);
rhs = T(141);
residual(103) = lhs - rhs;
lhs = y(104);
rhs = T(141);
residual(104) = lhs - rhs;
lhs = y(105);
rhs = T(141);
residual(105) = lhs - rhs;
lhs = y(106);
rhs = T(141);
residual(106) = lhs - rhs;
lhs = y(107);
rhs = T(141);
residual(107) = lhs - rhs;
lhs = y(108);
rhs = T(141);
residual(108) = lhs - rhs;
lhs = y(109);
rhs = T(141);
residual(109) = lhs - rhs;
lhs = y(110);
rhs = T(141);
residual(110) = lhs - rhs;
lhs = y(111);
rhs = T(141);
residual(111) = lhs - rhs;
lhs = y(112);
rhs = T(141);
residual(112) = lhs - rhs;
lhs = y(113);
rhs = T(141);
residual(113) = lhs - rhs;
lhs = y(114);
rhs = T(141);
residual(114) = lhs - rhs;
lhs = y(115);
rhs = T(141);
residual(115) = lhs - rhs;
lhs = y(116);
rhs = T(141);
residual(116) = lhs - rhs;
lhs = y(117);
rhs = T(141);
residual(117) = lhs - rhs;
lhs = y(118);
rhs = T(141);
residual(118) = lhs - rhs;
lhs = y(119);
rhs = T(141);
residual(119) = lhs - rhs;
lhs = y(120);
rhs = T(141);
residual(120) = lhs - rhs;
lhs = y(121);
rhs = T(141);
residual(121) = lhs - rhs;
lhs = y(122);
rhs = T(141);
residual(122) = lhs - rhs;
lhs = y(123);
rhs = T(141);
residual(123) = lhs - rhs;
lhs = y(124);
rhs = T(141);
residual(124) = lhs - rhs;
lhs = y(125);
rhs = T(141);
residual(125) = lhs - rhs;
lhs = y(126);
rhs = T(141);
residual(126) = lhs - rhs;
lhs = y(127);
rhs = T(141);
residual(127) = lhs - rhs;
lhs = y(128);
rhs = T(141);
residual(128) = lhs - rhs;
lhs = y(129);
rhs = T(141);
residual(129) = lhs - rhs;
lhs = y(130);
rhs = T(141);
residual(130) = lhs - rhs;
lhs = y(131);
rhs = T(141);
residual(131) = lhs - rhs;
lhs = y(132);
rhs = T(141);
residual(132) = lhs - rhs;
lhs = y(133);
rhs = T(141);
residual(133) = lhs - rhs;
lhs = y(134);
rhs = T(141);
residual(134) = lhs - rhs;
lhs = y(135);
rhs = y(47);
residual(135) = lhs - rhs;
lhs = y(136);
rhs = y(47);
residual(136) = lhs - rhs;
lhs = y(137);
rhs = y(47);
residual(137) = lhs - rhs;
lhs = y(138);
rhs = y(47);
residual(138) = lhs - rhs;
lhs = y(139);
rhs = y(47);
residual(139) = lhs - rhs;
lhs = y(140);
rhs = y(47);
residual(140) = lhs - rhs;
lhs = y(141);
rhs = y(47);
residual(141) = lhs - rhs;
lhs = y(142);
rhs = y(46);
residual(142) = lhs - rhs;
lhs = y(143);
rhs = y(46);
residual(143) = lhs - rhs;
lhs = y(144);
rhs = y(46);
residual(144) = lhs - rhs;
lhs = y(145);
rhs = y(46);
residual(145) = lhs - rhs;
lhs = y(146);
rhs = y(46);
residual(146) = lhs - rhs;
lhs = y(147);
rhs = y(46);
residual(147) = lhs - rhs;
lhs = y(148);
rhs = y(45);
residual(148) = lhs - rhs;
lhs = y(149);
rhs = y(45);
residual(149) = lhs - rhs;
lhs = y(150);
rhs = y(45);
residual(150) = lhs - rhs;
lhs = y(151);
rhs = y(45);
residual(151) = lhs - rhs;
lhs = y(152);
rhs = y(45);
residual(152) = lhs - rhs;
lhs = y(153);
rhs = y(44);
residual(153) = lhs - rhs;
lhs = y(154);
rhs = y(44);
residual(154) = lhs - rhs;
lhs = y(155);
rhs = y(44);
residual(155) = lhs - rhs;
lhs = y(156);
rhs = y(44);
residual(156) = lhs - rhs;
lhs = y(157);
rhs = y(43);
residual(157) = lhs - rhs;
lhs = y(158);
rhs = y(43);
residual(158) = lhs - rhs;
lhs = y(159);
rhs = y(43);
residual(159) = lhs - rhs;
lhs = y(160);
rhs = y(42);
residual(160) = lhs - rhs;
lhs = y(161);
rhs = y(42);
residual(161) = lhs - rhs;
lhs = y(162);
rhs = y(41);
residual(162) = lhs - rhs;
lhs = y(163);
rhs = y(36);
residual(163) = lhs - rhs;
lhs = y(164);
rhs = y(36);
residual(164) = lhs - rhs;
lhs = y(165);
rhs = y(36);
residual(165) = lhs - rhs;
lhs = y(166);
rhs = y(36);
residual(166) = lhs - rhs;
lhs = y(167);
rhs = y(36);
residual(167) = lhs - rhs;
lhs = y(168);
rhs = y(36);
residual(168) = lhs - rhs;
lhs = y(169);
rhs = y(36);
residual(169) = lhs - rhs;
lhs = y(170);
rhs = y(36);
residual(170) = lhs - rhs;
lhs = y(171);
rhs = y(36);
residual(171) = lhs - rhs;
lhs = y(172);
rhs = y(36);
residual(172) = lhs - rhs;
lhs = y(173);
rhs = y(36);
residual(173) = lhs - rhs;
lhs = y(174);
rhs = y(36);
residual(174) = lhs - rhs;
lhs = y(175);
rhs = y(36);
residual(175) = lhs - rhs;
lhs = y(176);
rhs = y(36);
residual(176) = lhs - rhs;
lhs = y(177);
rhs = y(36);
residual(177) = lhs - rhs;
lhs = y(178);
rhs = y(36);
residual(178) = lhs - rhs;
lhs = y(179);
rhs = y(36);
residual(179) = lhs - rhs;
lhs = y(180);
rhs = y(36);
residual(180) = lhs - rhs;
lhs = y(181);
rhs = y(36);
residual(181) = lhs - rhs;
lhs = y(182);
rhs = y(36);
residual(182) = lhs - rhs;
lhs = y(183);
rhs = y(36);
residual(183) = lhs - rhs;
lhs = y(184);
rhs = y(36);
residual(184) = lhs - rhs;
lhs = y(185);
rhs = y(36);
residual(185) = lhs - rhs;
lhs = y(186);
rhs = y(36);
residual(186) = lhs - rhs;
lhs = y(187);
rhs = y(36);
residual(187) = lhs - rhs;
lhs = y(188);
rhs = y(36);
residual(188) = lhs - rhs;
lhs = y(189);
rhs = y(36);
residual(189) = lhs - rhs;
lhs = y(190);
rhs = y(36);
residual(190) = lhs - rhs;
lhs = y(191);
rhs = y(36);
residual(191) = lhs - rhs;
lhs = y(192);
rhs = y(36);
residual(192) = lhs - rhs;
lhs = y(193);
rhs = y(36);
residual(193) = lhs - rhs;
lhs = y(194);
rhs = y(36);
residual(194) = lhs - rhs;
lhs = y(195);
rhs = y(36);
residual(195) = lhs - rhs;
lhs = y(196);
rhs = y(36);
residual(196) = lhs - rhs;
lhs = y(197);
rhs = y(36);
residual(197) = lhs - rhs;
lhs = y(198);
rhs = y(36);
residual(198) = lhs - rhs;
lhs = y(199);
rhs = y(36);
residual(199) = lhs - rhs;
lhs = y(200);
rhs = y(36);
residual(200) = lhs - rhs;
lhs = y(201);
rhs = y(49);
residual(201) = lhs - rhs;
lhs = y(202);
rhs = y(49);
residual(202) = lhs - rhs;
lhs = y(203);
rhs = y(49);
residual(203) = lhs - rhs;
lhs = y(204);
rhs = y(49);
residual(204) = lhs - rhs;
lhs = y(205);
rhs = y(49);
residual(205) = lhs - rhs;
lhs = y(206);
rhs = y(49);
residual(206) = lhs - rhs;
lhs = y(207);
rhs = y(49);
residual(207) = lhs - rhs;
lhs = y(208);
rhs = y(49);
residual(208) = lhs - rhs;
lhs = y(209);
rhs = y(49);
residual(209) = lhs - rhs;
lhs = y(210);
rhs = y(49);
residual(210) = lhs - rhs;
lhs = y(211);
rhs = y(49);
residual(211) = lhs - rhs;
lhs = y(212);
rhs = y(49);
residual(212) = lhs - rhs;
lhs = y(213);
rhs = y(49);
residual(213) = lhs - rhs;
lhs = y(214);
rhs = y(49);
residual(214) = lhs - rhs;
lhs = y(215);
rhs = y(49);
residual(215) = lhs - rhs;
lhs = y(216);
rhs = y(49);
residual(216) = lhs - rhs;
lhs = y(217);
rhs = y(49);
residual(217) = lhs - rhs;
lhs = y(218);
rhs = y(49);
residual(218) = lhs - rhs;
lhs = y(219);
rhs = y(49);
residual(219) = lhs - rhs;
lhs = y(220);
rhs = y(49);
residual(220) = lhs - rhs;
lhs = y(221);
rhs = y(49);
residual(221) = lhs - rhs;
lhs = y(222);
rhs = y(49);
residual(222) = lhs - rhs;
lhs = y(223);
rhs = y(49);
residual(223) = lhs - rhs;
lhs = y(224);
rhs = y(49);
residual(224) = lhs - rhs;
lhs = y(225);
rhs = y(49);
residual(225) = lhs - rhs;
lhs = y(226);
rhs = y(49);
residual(226) = lhs - rhs;
lhs = y(227);
rhs = y(49);
residual(227) = lhs - rhs;
lhs = y(228);
rhs = y(49);
residual(228) = lhs - rhs;
lhs = y(229);
rhs = y(49);
residual(229) = lhs - rhs;
lhs = y(230);
rhs = y(49);
residual(230) = lhs - rhs;
lhs = y(231);
rhs = y(49);
residual(231) = lhs - rhs;
lhs = y(232);
rhs = y(49);
residual(232) = lhs - rhs;
lhs = y(233);
rhs = y(49);
residual(233) = lhs - rhs;
lhs = y(234);
rhs = y(49);
residual(234) = lhs - rhs;
lhs = y(235);
rhs = y(49);
residual(235) = lhs - rhs;
lhs = y(236);
rhs = y(49);
residual(236) = lhs - rhs;
lhs = y(237);
rhs = y(49);
residual(237) = lhs - rhs;
lhs = y(238);
rhs = y(24);
residual(238) = lhs - rhs;
lhs = y(239);
rhs = y(24);
residual(239) = lhs - rhs;
lhs = y(240);
rhs = y(24);
residual(240) = lhs - rhs;
lhs = y(241);
rhs = y(24);
residual(241) = lhs - rhs;
lhs = y(242);
rhs = y(24);
residual(242) = lhs - rhs;
lhs = y(243);
rhs = y(24);
residual(243) = lhs - rhs;
lhs = y(244);
rhs = y(24);
residual(244) = lhs - rhs;
lhs = y(245);
rhs = y(24);
residual(245) = lhs - rhs;
lhs = y(246);
rhs = y(24);
residual(246) = lhs - rhs;
lhs = y(247);
rhs = y(24);
residual(247) = lhs - rhs;
lhs = y(248);
rhs = y(24);
residual(248) = lhs - rhs;
lhs = y(249);
rhs = y(24);
residual(249) = lhs - rhs;
lhs = y(250);
rhs = y(24);
residual(250) = lhs - rhs;
lhs = y(251);
rhs = y(24);
residual(251) = lhs - rhs;
lhs = y(252);
rhs = y(24);
residual(252) = lhs - rhs;
lhs = y(253);
rhs = y(24);
residual(253) = lhs - rhs;
lhs = y(254);
rhs = y(24);
residual(254) = lhs - rhs;
lhs = y(255);
rhs = y(24);
residual(255) = lhs - rhs;
lhs = y(256);
rhs = y(24);
residual(256) = lhs - rhs;
lhs = y(257);
rhs = y(24);
residual(257) = lhs - rhs;
lhs = y(258);
rhs = y(24);
residual(258) = lhs - rhs;
lhs = y(259);
rhs = y(24);
residual(259) = lhs - rhs;
lhs = y(260);
rhs = y(24);
residual(260) = lhs - rhs;
lhs = y(261);
rhs = y(24);
residual(261) = lhs - rhs;
lhs = y(262);
rhs = y(24);
residual(262) = lhs - rhs;
lhs = y(263);
rhs = y(24);
residual(263) = lhs - rhs;
lhs = y(264);
rhs = y(24);
residual(264) = lhs - rhs;
lhs = y(265);
rhs = y(24);
residual(265) = lhs - rhs;
lhs = y(266);
rhs = y(24);
residual(266) = lhs - rhs;
lhs = y(267);
rhs = y(24);
residual(267) = lhs - rhs;
lhs = y(268);
rhs = y(24);
residual(268) = lhs - rhs;
lhs = y(269);
rhs = y(24);
residual(269) = lhs - rhs;
lhs = y(270);
rhs = y(24);
residual(270) = lhs - rhs;
lhs = y(271);
rhs = y(24);
residual(271) = lhs - rhs;
lhs = y(272);
rhs = y(24);
residual(272) = lhs - rhs;
lhs = y(273);
rhs = y(24);
residual(273) = lhs - rhs;
lhs = y(274);
rhs = y(24);
residual(274) = lhs - rhs;
lhs = y(275);
rhs = y(19);
residual(275) = lhs - rhs;
lhs = y(276);
rhs = y(19);
residual(276) = lhs - rhs;
lhs = y(277);
rhs = y(19);
residual(277) = lhs - rhs;
lhs = y(278);
rhs = y(19);
residual(278) = lhs - rhs;
lhs = y(279);
rhs = y(19);
residual(279) = lhs - rhs;
lhs = y(280);
rhs = y(19);
residual(280) = lhs - rhs;
lhs = y(281);
rhs = y(19);
residual(281) = lhs - rhs;
lhs = y(282);
rhs = y(19);
residual(282) = lhs - rhs;
lhs = y(283);
rhs = y(19);
residual(283) = lhs - rhs;
lhs = y(284);
rhs = y(19);
residual(284) = lhs - rhs;
lhs = y(285);
rhs = y(19);
residual(285) = lhs - rhs;
lhs = y(286);
rhs = y(19);
residual(286) = lhs - rhs;
lhs = y(287);
rhs = y(19);
residual(287) = lhs - rhs;
lhs = y(288);
rhs = y(19);
residual(288) = lhs - rhs;
lhs = y(289);
rhs = y(19);
residual(289) = lhs - rhs;
lhs = y(290);
rhs = y(19);
residual(290) = lhs - rhs;
lhs = y(291);
rhs = y(19);
residual(291) = lhs - rhs;
lhs = y(292);
rhs = y(19);
residual(292) = lhs - rhs;
lhs = y(293);
rhs = y(19);
residual(293) = lhs - rhs;
lhs = y(294);
rhs = y(19);
residual(294) = lhs - rhs;
lhs = y(295);
rhs = y(19);
residual(295) = lhs - rhs;
lhs = y(296);
rhs = y(19);
residual(296) = lhs - rhs;
lhs = y(297);
rhs = y(19);
residual(297) = lhs - rhs;
lhs = y(298);
rhs = y(19);
residual(298) = lhs - rhs;
lhs = y(299);
rhs = y(19);
residual(299) = lhs - rhs;
lhs = y(300);
rhs = y(19);
residual(300) = lhs - rhs;
lhs = y(301);
rhs = y(19);
residual(301) = lhs - rhs;
lhs = y(302);
rhs = y(19);
residual(302) = lhs - rhs;
lhs = y(303);
rhs = y(19);
residual(303) = lhs - rhs;
lhs = y(304);
rhs = y(19);
residual(304) = lhs - rhs;
lhs = y(305);
rhs = y(19);
residual(305) = lhs - rhs;
lhs = y(306);
rhs = y(19);
residual(306) = lhs - rhs;
lhs = y(307);
rhs = y(19);
residual(307) = lhs - rhs;
lhs = y(308);
rhs = y(19);
residual(308) = lhs - rhs;
lhs = y(309);
rhs = y(19);
residual(309) = lhs - rhs;
lhs = y(310);
rhs = y(19);
residual(310) = lhs - rhs;
lhs = y(311);
rhs = y(19);
residual(311) = lhs - rhs;
lhs = y(312);
rhs = y(33);
residual(312) = lhs - rhs;
lhs = y(313);
rhs = y(33);
residual(313) = lhs - rhs;
lhs = y(314);
rhs = y(33);
residual(314) = lhs - rhs;
lhs = y(315);
rhs = y(33);
residual(315) = lhs - rhs;
lhs = y(316);
rhs = y(33);
residual(316) = lhs - rhs;
lhs = y(317);
rhs = y(33);
residual(317) = lhs - rhs;
lhs = y(318);
rhs = y(33);
residual(318) = lhs - rhs;
lhs = y(319);
rhs = y(33);
residual(319) = lhs - rhs;
lhs = y(320);
rhs = y(33);
residual(320) = lhs - rhs;
lhs = y(321);
rhs = y(33);
residual(321) = lhs - rhs;
lhs = y(322);
rhs = y(33);
residual(322) = lhs - rhs;
lhs = y(323);
rhs = y(33);
residual(323) = lhs - rhs;
lhs = y(324);
rhs = y(33);
residual(324) = lhs - rhs;
lhs = y(325);
rhs = y(33);
residual(325) = lhs - rhs;
lhs = y(326);
rhs = y(33);
residual(326) = lhs - rhs;
lhs = y(327);
rhs = y(33);
residual(327) = lhs - rhs;
lhs = y(328);
rhs = y(33);
residual(328) = lhs - rhs;
lhs = y(329);
rhs = y(33);
residual(329) = lhs - rhs;
lhs = y(330);
rhs = y(33);
residual(330) = lhs - rhs;
lhs = y(331);
rhs = y(33);
residual(331) = lhs - rhs;
lhs = y(332);
rhs = y(33);
residual(332) = lhs - rhs;
lhs = y(333);
rhs = y(33);
residual(333) = lhs - rhs;
lhs = y(334);
rhs = y(33);
residual(334) = lhs - rhs;
lhs = y(335);
rhs = y(33);
residual(335) = lhs - rhs;
lhs = y(336);
rhs = y(33);
residual(336) = lhs - rhs;
lhs = y(337);
rhs = y(33);
residual(337) = lhs - rhs;
lhs = y(338);
rhs = y(33);
residual(338) = lhs - rhs;
lhs = y(339);
rhs = y(33);
residual(339) = lhs - rhs;
lhs = y(340);
rhs = y(33);
residual(340) = lhs - rhs;
lhs = y(341);
rhs = y(33);
residual(341) = lhs - rhs;
lhs = y(342);
rhs = y(33);
residual(342) = lhs - rhs;
lhs = y(343);
rhs = y(33);
residual(343) = lhs - rhs;
lhs = y(344);
rhs = y(33);
residual(344) = lhs - rhs;
lhs = y(345);
rhs = y(33);
residual(345) = lhs - rhs;
lhs = y(346);
rhs = y(33);
residual(346) = lhs - rhs;
lhs = y(347);
rhs = y(33);
residual(347) = lhs - rhs;
lhs = y(348);
rhs = y(33);
residual(348) = lhs - rhs;
lhs = y(349);
rhs = y(33);
residual(349) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end