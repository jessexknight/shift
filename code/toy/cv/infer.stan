data {
  int N;  // num individuals
  int Ne; // num exposures
  int Nk; // num strata = 2^Ne
  real Py; // prevalence of outcome
  vector[Ne] Pe; // prevalence of exposures
  vector[Ne] OR; // odds ratios: outcome vs exposure
  array[Ne,Nk/2] int e0; // exposure indices: e=0
  array[Ne,Nk/2] int e1; // exposure indices: e=1
}

parameters {
  simplex[Nk*2] p; // latent probability: exposure + outcome
}

transformed parameters {
  // expected counts p[e*y] -> kye[y][e]
  array[2] vector[Nk] kye;
  kye[1] = N * p[0*Nk+1:Nk*1]; // 1st half
  kye[2] = N * p[1*Nk+1:Nk*2]; // 2nd half
  // total counts per (strata,outcome)
  array[Ne,4] real k4;
  for (i in 1:Ne){
    k4[i,1] = sum(kye[1,e0[i]]); // e=0,y=0
    k4[i,2] = sum(kye[1,e1[i]]); // e=1,y=0
    k4[i,3] = sum(kye[2,e0[i]]); // e=0,y=1
    k4[i,4] = sum(kye[2,e1[i]]); // e=1,y=1
  }
}

model {
  vector[Ne] LORmu; // expected log(OR) mean
  vector[Ne] LORse; // expected log(OR) SD
  // observed outcome prevalence
  Py ~ beta(k4[1,3]+k4[1,4], k4[1,1]+k4[1,2]);
  for (i in 1:Ne){
    // observed exposure prevalence
    Pe[i] ~ beta(k4[i,2]+k4[i,4], k4[i,1]+k4[i,3]);
    // observed odds ratios
    LORmu[i] = log((k4[i,1]*k4[i,4]) / (k4[i,2]*k4[i,3]));
    LORse[i] = sqrt(1/k4[i,1]+1/k4[i,2]+1/k4[i,3]+1/k4[i,4]);
    log(OR[i]) ~ normal(LORmu[i],LORse[i]);
  }
}

generated quantities {
  vector[Nk] ke  = kye[1] + kye[2];      // counts per exposure
  vector[Nk] pye = kye[2] ./ ke;         // outcome prev per exposure
  real pym = sum(ke .* pye) / N;         // weighted mean prevalence
  real pyv = sum(ke .* (pye-pym)^2) / N; // weighted var prevalence
  real pys = sqrt(pyv);                  // weighted SD prevalence
  real cv = pys / pym;                   // coefficient of variation
}
