char *RealName(int idhep) {
  // Got it from http://root.cern.ch/root/html/src/TAttParticle.cxx.html
  switch (idhep) {
  case  0 : return "frag.";
  case  1 : return "down";
  case  -1 : return "down bar";
  case  2 : return "up";
  case  -2 : return "up bar";
  case  3 : return "strange";
  case  -3 : return "strange bar";
  case  4 : return "charm";
  case  -4 : return "charm bar";
  case  5 : return "bottom";
  case  -5 : return "bottom bar";
  case  6 : return "top";
  case  -6 : return "top bar";
  case  21 : return "gluon";
  case  7 : return "Searches0";
  case  11 : return "e-";
  case  -11 : return "e+";
  case  12 : return "nu(e)";
  case  -12 : return "nu(e) bar";
  case  13 : return "mu-";
  case  -13 : return "mu+";
  case  14 : return "nu(mu)";
  case  -14 : return "nu(mu) bar";
  case  15 : return "tau-";
  case  -15 : return "tau+";
  case  16 : return "nu(tau)";
  case  -16 : return "nu(tau) bar";
  case  22 : return "gamma";
  case  23 : return "Z0";
  case  24 : return "W+";
  case  -24 : return "W-";
  case  111 : return "pi0";
  case  113 : return "rho(770)0";
  case  115 : return "a(2)(1320)0";
  case  117 : return "rho(3)(1690)0";
  case  130 : return "K(L)0";
  case  211 : return "pi+";
  case  -211 : return "pi-";
  case  213 : return "rho(770)+";
  case  -213 : return "rho(770)-";
  case  215 : return "a(2)(1320)+";
  case  -215 : return "a(2)(1320)-";
  case  217 : return "rho(3)(1690)+";
  case  -217 : return "rho(3)(1690)-";
  case  221 : return "eta0";
  case  223 : return "omega(782)0";
  case  225 : return "f(2)(1270)0";
  case  227 : return "omega(3)(1670)0";
  case  229 : return "f(4)(2050)0";
  case  310 : return "K(S)0";
  case  311 : return "K0";
  case  -311 : return "K0 bar";
  case  313 : return "K*(892)0";
  case  -313 : return "K*(892)0 bar";
  case  315 : return "K(2)*(1430)0";
  case  -315 : return "K(2)*(1430)0 bar";
  case  317 : return "K(3)*(1780)0";
  case  -317 : return "K(3)*(1780)0 bar";
  case  319 : return "K(4)*(2045)0";
  case  -319 : return "K(4)*(2045)0 bar";
  case  321 : return "K+";
  case  -321 : return "K-";
  case  323 : return "K*(892)+";
  case  -323 : return "K*(892)-";
  case  325 : return "K(2)*(1430)+";
  case  -325 : return "K(2)*(1430)-";
  case  327 : return "K(3)*(1780)+";
  case  -327 : return "K(3)*(1780)-";
  case  329 : return "K(4)*(2045)+";
  case  -329 : return "K(4)*(2045)-";
  case  331 : return "eta'(958)0";
  case  333 : return "phi(1020)0";
  case  335 : return "f(2)'(1525)0";
  case  337 : return "phi(3)(1850)0";
  case  411 : return "D+";
  case  -411 : return "D-";
  case  413 : return "D*(2010)+";
  case  -413 : return "D*(2010)-";
  case  415 : return "D(2)*(2460)+";
  case  -415 : return "D(2)*(2460)-";
  case  421 : return "D0";
  case  -421 : return "D0~";
  case  423 : return "D*(2007)0";
  case  -423 : return "D*(2007)0~";
  case  425 : return "D(2)*(2460)0";
  case  431 : return "D(s)+";
  case  -431 : return "D(s)-";
  case  433 : return "D(s)*+";
  case  -433 : return "D(s)*-";
  case  441 : return "eta(c)(1S)0";
  case  443 : return "J/psi(1S)0";
  case  445 : return "chi(c2)(1P)0";
  case  511 : return "B0";
  case  -511 : return "B0~";
  case  513 : return "B*0";
  case  521 : return "B+";
  case  -521 : return "B-";
  case  523 : return "B*+";
  case  -523 : return "B*-";
  case  531 : return "B(s)0";
  case  551 : return "chi(b0)(1P)0";
  case  553 : return "Upsilon(1S)0";
  case  555 : return "chi(b2)(1P)0";
  case  1112 : return "Delta(1620)-";
  case  -1112 : return "Delta(1620)+ bar";
  case  1114 : return "Delta(1232)-";
  case  -1114 : return "Delta(1232)+ bar";
  case  1116 : return "Delta(1905)-";
  case  -1116 : return "Delta(1905)+ bar";
  case  1118 : return "Delta(1950)-";
  case  -1118 : return "Delta(1950)+ bar";
  case  1212 : return "Delta(1620)0";
  case  -1212 : return "Delta(1620)0 bar";
  case  1214 : return "N(1520)0";
  case  -1214 : return "N(1520)0 bar";
  case  1216 : return "Delta(1905)0";
  case  -1216 : return "Delta(1905)0 bar";
  case  1218 : return "N(2190)0";
  case  -1218 : return "N(2190)0 bar";
  case  2112 : return "n";
  case  -2112 : return "n bar";
  case  2114 : return "Delta(1232)0";
  case  -2114 : return "Delta(1232)0 bar";
  case  2116 : return "N(1675)0";
  case  -2116 : return "N(1675)0 bar";
  case  2118 : return "Delta(1950)0";
  case  -2118 : return "Delta(1950)0 bar";
  case  2122 : return "Delta(1620)+";
  case  -2122 : return "Delta(1620)- bar";
  case  2124 : return "N(1520)+";
  case  -2124 : return "N(1520)- bar";
  case  2126 : return "Delta(1905)+";
  case  -2126 : return "Delta(1905)- bar";
  case  2128 : return "N(2190)+";
  case  -2128 : return "N(2190)- bar";
  case  2212 : return "p";
  case  -2212 : return "p bar";
  case  2214 : return "Delta(1232)+";
  case  -2214 : return "Delta(1232)- bar";
  case  2216 : return "N(1675)+";
  case  -2216 : return "N(1675)- bar";
  case  2218 : return "Delta(1950)+";
  case  -2218 : return "Delta(1950)- bar";
  case  2222 : return "Delta(1620)++";
  case  -2222 : return "Delta(1620)-- bar";
  case  2224 : return "Delta(1232)++";
  case  -2224 : return "Delta(1232)-- bar";
  case  2226 : return "Delta(1905)++";
  case  -2226 : return "Delta(1905)-- bar";
  case  2228 : return "Delta(1950)++";
  case  -2228 : return "Delta(1950)-- bar";
  case  3112 : return "Sigma-";
  case  -3112 : return "Sigma+ bar";
  case  3114 : return "Sigma(1385)-";
  case  -3114 : return "Sigma(1385)+ bar";
  case  3116 : return "Sigma(1775)-";
  case  -3116 : return "Sigma(1775)+ bar";
  case  3118 : return "Sigma(2030)-";
  case  -3118 : return "Sigma(2030)+ bar";
  case  3122 : return "Lambda0";
  case  -3122 : return "Lambda0 bar";
  case  3124 : return "Lambda(1520)0";
  case  -3124 : return "Lambda(1520)0 bar";
  case  3126 : return "Lambda(1820)0";
  case  -3126 : return "Lambda(1820)0 bar";
  case  3128 : return "Lambda(2100)0";
  case  -3128 : return "Lambda(2100)0 bar";
  case  3212 : return "Sigma0";
  case  -3212 : return "Sigma0 bar";
  case  3214 : return "Sigma(1385)0";
  case  -3214 : return "Sigma(1385)0 bar";
  case  3216 : return "Sigma(1775)0";
  case  -3216 : return "Sigma(1775)0 bar";
  case  3218 : return "Sigma(2030)0";
  case  -3218 : return "Sigma(2030)0 bar";
  case  3222 : return "Sigma+";
  case  -3222 : return "Sigma- bar";
  case  3224 : return "Sigma(1385)+";
  case  -3224 : return "Sigma(1385)- bar";
  case  3226 : return "Sigma(1775)+";
  case  -3226 : return "Sigma(1775)- bar";
  case  3228 : return "Sigma(2030)+";
  case  -3228 : return "Sigma(2030)- bar";
  case  3312 : return "Xi-";
  case  -3312 : return "Xi+ bar";
  case  3314 : return "Xi(1530)-";
  case  -3314 : return "Xi(1530)+ bar";
  case  3322 : return "Xi0";
  case  -3322 : return "Xi0 bar";
  case  3324 : return "Xi(1530)0";
  case  -3324 : return "Xi(1530)0 bar";
  case  3334 : return "Omega-";
  case  -3334 : return "Omega+ bar";
  case  4112 : return "Sigma(c)(2455)0";
  case  -4112 : return "Sigma(c)(2455)0 bar";
  case  4122 : return "Lambda(c)+";
  case  -4122 : return "Lambda(c)- bar";
  case  4212 : return "Sigma(c)(2455)+";
  case  -4212 : return "Sigma(c)(2455)- bar";
  case  4222 : return "Sigma(c)(2455)++";
  case  -4222 : return "Sigma(c)(2455)-- bar";
  case  4312 : return "Xi(c)0";
  case  -4312 : return "Xi(c)0 bar";
  case  4322 : return "Xi(c)+";
  case  -4322 : return "Xi(c)- bar";
  case  5122 : return "Lambda(b)0";
  case  -5122 : return "Lambda(b)0 bar";
  case  10111 : return "a(0)(980)0";
  case  10113 : return "b(1)(1235)0";
  case  10115 : return "pi(2)(1670)0";
  case  10211 : return "a(0)(980)+";
  case  -10211 : return "a(0)(980)-";
  case  10213 : return "b(1)(1235)+";
  case  -10213 : return "b(1)(1235)-";
  case  10215 : return "pi(2)(1670)+";
  case  -10215 : return "pi(2)(1670)-";
  case  10221 : return "f(0)(980)0";
  case  10223 : return "h(1)(1170)0";
  case  10311 : return "K(0)*(1430)0";
  case  -10311 : return "K(0)*(1430)0 bar";
  case  10313 : return "K(1)(1270)0";
  case  -10313 : return "K(1)(1270)0 bar";
  case  10315 : return "K(2)(1770)0";
  case  -10315 : return "K(2)(1770)0 bar";
  case  10321 : return "K(0)*(1430)+";
  case  -10321 : return "K(0)*(1430)-";
  case  10323 : return "K(1)(1270)+";
  case  -10323 : return "K(1)(1270)-";
  case  10325 : return "K(2)(1770)+";
  case  -10325 : return "K(2)(1770)-";
  case  10333 : return "phi(1680)0";
  case  10423 : return "D(1)(2420)0";
  case  10433 : return "D(s1)(2536)+";
  case  -10433 : return "D(s1)(2536)-";
  case  10441 : return "chi(c0)(1P)0";
  case  10443 : return "chi(c1)(1P)0";
  case  10551 : return "chi(b0)(2P)0";
  case  10553 : return "chi(b1)(1P)0";
  case  10555 : return "chi(b2)(2P)0";
  case  11112 : return "Delta(1900)-";
  case  -11112 : return "Delta(1900)+ bar";
  case  11114 : return "Delta(1700)-";
  case  -11114 : return "Delta(1700)+ bar";
  case  11116 : return "Delta(1930)-";
  case  -11116 : return "Delta(1930)+ bar";
  case  11212 : return "Delta(1900)0";
  case  -11212 : return "Delta(1900)0 bar";
  case  11216 : return "Delta(1930)0";
  case  -11216 : return "Delta(1930)0 bar";
  case  12112 : return "N(1440)0";
  case  -12112 : return "N(1440)0 bar";
  case  12114 : return "Delta(1700)0";
  case  -12114 : return "Delta(1700)0 bar";
  case  12116 : return "N(1680)0";
  case  -12116 : return "N(1680)0 bar";
  case  12122 : return "Delta(1900)+";
  case  -12122 : return "Delta(1900)- bar";
  case  12126 : return "Delta(1930)+";
  case  -12126 : return "Delta(1930)- bar";
  case  12212 : return "N(1440)+";
  case  -12212 : return "N(1440)- bar";
  case  12214 : return "Delta(1700)+";
  case  -12214 : return "Delta(1700)- bar";
  case  12216 : return "N(1680)+";
  case  -12216 : return "N(1680)- bar";
  case  12222 : return "Delta(1900)++";
  case  -12222 : return "Delta(1900)-- bar";
  case  12224 : return "Delta(1700)++";
  case  -12224 : return "Delta(1700)-- bar";
  case  12226 : return "Delta(1930)++";
  case  -12226 : return "Delta(1930)-- bar";
  case  13112 : return "Sigma(1660)-";
  case  -13112 : return "Sigma(1660)+ bar";
  case  13114 : return "Sigma(1670)-";
  case  -13114 : return "Sigma(1670)+ bar";
  case  13116 : return "Sigma(1915)-";
  case  -13116 : return "Sigma(1915)+ bar";
  case  13122 : return "Lambda(1405)0";
  case  -13122 : return "Lambda(1405)0 bar";
  case  13124 : return "Lambda(1690)0";
  case  -13124 : return "Lambda(1690)0 bar";
  case  13126 : return "Lambda(1830)0";
  case  -13126 : return "Lambda(1830)0 bar";
  case  13212 : return "Sigma(1660)0";
  case  -13212 : return "Sigma(1660)0 bar";
  case  13214 : return "Sigma(1670)0";
  case  -13214 : return "Sigma(1670)0 bar";
  case  13216 : return "Sigma(1915)0";
  case  -13216 : return "Sigma(1915)0 bar";
  case  13222 : return "Sigma(1660)+";
  case  -13222 : return "Sigma(1660)- bar";
  case  13224 : return "Sigma(1670)+";
  case  -13224 : return "Sigma(1670)- bar";
  case  13226 : return "Sigma(1915)+";
  case  -13226 : return "Sigma(1915)- bar";
  case  13314 : return "Xi(1820)-";
  case  -13314 : return "Xi(1820)+ bar";
  case  13324 : return "Xi(1820)0";
  case  -13324 : return "Xi(1820)0 bar";
  case  20111 : return "pi(1300)0";
  case  20113 : return "a(1)(1260)0";
  case  20211 : return "pi(1300)+";
  case  -20211 : return "pi(1300)-";
  case  20213 : return "a(1)(1260)+";
  case  -20213 : return "a(1)(1260)-";
  case  20221 : return "eta(1295)0";
  case  20223 : return "f(1)(1285)0";
  case  20225 : return "f(2)(2010)0";
  case  20313 : return "K(1)(1400)0";
  case  -20313 : return "K(1)(1400)0 bar";
  case  20315 : return "K(2)(1820)0";
  case  -20315 : return "K(2)(1820)0 bar";
  case  20323 : return "K(1)(1400)+";
  case  -20323 : return "K(1)(1400)-";
  case  20325 : return "K(2)(1820)+";
  case  -20325 : return "K(2)(1820)-";
  case  20443 : return "psi(2S)0";
  case  20553 : return "Upsilon(2S)0";
  case  21112 : return "Delta(1910)-";
  case  -21112 : return "Delta(1910)+ bar";
  case  21114 : return "Delta(1920)-";
  case  -21114 : return "Delta(1920)+ bar";
  case  21212 : return "Delta(1910)0";
  case  -21212 : return "Delta(1910)0 bar";
  case  21214 : return "N(1700)0";
  case  -21214 : return "N(1700)0 bar";
  case  22112 : return "N(1535)0";
  case  -22112 : return "N(1535)0 bar";
  case  22114 : return "Delta(1920)0";
  case  -22114 : return "Delta(1920)0 bar";
  case  22122 : return "Delta(1910)+";
  case  -22122 : return "Delta(1910)- bar";
  case  22124 : return "N(1700)+";
  case  -22124 : return "N(1700)- bar";
  case  22212 : return "N(1535)+";
  case  -22212 : return "N(1535)- bar";
  case  22214 : return "Delta(1920)+";
  case  -22214 : return "Delta(1920)- bar";
  case  22222 : return "Delta(1910)++";
  case  -22222 : return "Delta(1910)-- bar";
  case  22224 : return "Delta(1920)++";
  case  -22224 : return "Delta(1920)-- bar";
  case  23112 : return "Sigma(1750)-";
  case  -23112 : return "Sigma(1750)+ bar";
  case  23114 : return "Sigma(1940)-";
  case  -23114 : return "Sigma(1940)+ bar";
  case  23122 : return "Lambda(1600)0";
  case  -23122 : return "Lambda(1600)0 bar";
  case  23124 : return "Lambda(1890)0";
  case  -23124 : return "Lambda(1890)0 bar";
  case  23126 : return "Lambda(2110)0";
  case  -23126 : return "Lambda(2110)0 bar";
  case  23212 : return "Sigma(1750)0";
  case  -23212 : return "Sigma(1750)0 bar";
  case  23214 : return "Sigma(1940)0";
  case  -23214 : return "Sigma(1940)0 bar";
  case  23222 : return "Sigma(1750)+";
  case  -23222 : return "Sigma(1750)- bar";
  case  23224 : return "Sigma(1940)+";
  case  -23224 : return "Sigma(1940)- bar";
  case  30113 : return "rho(1700)0";
  case  30213 : return "rho(1700)+";
  case  -30213 : return "rho(1700)-";
  case  30223 : return "f(1)(1420)0";
  case  30225 : return "f(2)(2300)0";
  case  30313 : return "K*(1410)0";
  case  -30313 : return "K*(1410)0 bar";
  case  30323 : return "K*(1410)+";
  case  -30323 : return "K*(1410)-";
  case  30443 : return "psi(3770)0";
  case  30553 : return "Upsilon(3S)0";
  case  31114 : return "Delta(1600)-";
  case  -31114 : return "Delta(1600)+ bar";
  case  31214 : return "N(1720)0";
  case  -31214 : return "N(1720)0 bar";
  case  32112 : return "N(1650)0";
  case  -32112 : return "N(1650)0 bar";
  case  32114 : return "Delta(1600)0";
  case  -32114 : return "Delta(1600)0 bar";
  case  32124 : return "N(1720)+";
  case  -32124 : return "N(1720)- bar";
  case  32212 : return "N(1650)+";
  case  -32212 : return "N(1650)- bar";
  case  32214 : return "Delta(1600)+";
  case  -32214 : return "Delta(1600)- bar";
  case  32224 : return "Delta(1600)++";
  case  -32224 : return "Delta(1600)-- bar";
  case  33122 : return "Lambda(1670)0";
  case  -33122 : return "Lambda(1670)0 bar";
  case  40113 : return "rho(1450)0";
  case  40213 : return "rho(1450)+";
  case  -40213 : return "rho(1450)-";
  case  40221 : return "eta(1440)0";
  case  40223 : return "f(1)(1510)0";
  case  40225 : return "f(2)(2340)0";
  case  40313 : return "K*(1680)0";
  case  -40313 : return "K*(1680)0 bar";
  case  40323 : return "K*(1680)+";
  case  -40323 : return "K*(1680)-";
  case  40443 : return "psi(4040)0";
  case  40553 : return "Upsilon(4S)0";
  case  42112 : return "N(1710)0";
  case  -42112 : return "N(1710)0 bar";
  case  42212 : return "N(1710)+";
  case  -42212 : return "N(1710)- bar";
  case  43122 : return "Lambda(1800)0";
  case  -43122 : return "Lambda(1800)0 bar";
  case  50221 : return "f(0)(1590)0";
  case  50223 : return "omega(1420)0";
  case  50443 : return "psi(4160)0";
  case  50553 : return "Upsilon(10860)0";
  case  53122 : return "Lambda(1810)0";
  case  -53122 : return "Lambda(1810)0 bar";
  case  60221 : return "f(J)(1710)0";
  case  60223 : return "omega(1600)0";
  case  60443 : return "psi(4415)0";
  case  60553 : return "Upsilon(11020)0";
  case  70553 : return "chi(b1)(2P)0";
  case     25 : return "h^0";
  case     35 : return "H^0";
  case     36 : return "A^0";
  case     37 : return "H^+-";
  case     1000001 : return "tilde{d}_L";
  case     1000002 : return "tilde{u}_L";
  case     1000003 : return "tilde{s}_L";
  case     1000004 : return "tilde{c}_L";
  case     1000005 : return "tilde{b}_1";
  case     1000006 : return "tilde{t}_1";
  case     2000001 : return "tilde{d}_R";
  case     2000002 : return "tilde{u}_R";
  case     2000003 : return "tilde{s}_R";
  case     2000004 : return "tilde{c}_R";
  case     2000005 : return "tilde{b}_2";
  case     2000006 : return "tilde{t}_2";
  case     1000011 : return "tilde{e}_L";
  case     1000012 : return "tilde{nu}_L";
  case     1000013 : return "tilde{mu}_L";
  case     1000015 : return "tilde{tau}_1";
  case     2000011 : return "tilde{e}_R";
  case     2000012 : return "tilde{nu}_R";
  case     2000013 : return "tilde{mu}_R";
  case     2000015 : return "tilde{tau}_2";
  case     1000021 : return "tilde{g}";
  case     1000022 : return "tilde{chi}_1^0";
  case     1000023 : return "tilde{chi}_2^0";
  case     1000024 : return "tilde{chi}_1^+-";
  case     1000025 : return "tilde{chi}_3^0";
  case     1000035 : return "tilde{chi}_4^0";
  case     1000037 : return "tilde{chi}_2^+-";


  default : return "Unknown";      // isajet or pdg number does not exist
 }
}
