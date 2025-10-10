data {
  int N;  // num individuals
  int Ne; // num exposures
  int Ns; // num strata = 2^Ne
  real Py; // prevalence of outcome
  vector[Ne] Pe; // prevalence of exposures
  vector[Ne] OR; // odds ratios: outcome vs exposure
  array[Ne,Ns/2] int se0; // exposure indices: e=0
  array[Ne,Ns/2] int se1; // exposure indices: e=1
}

parameters {
  simplex[Ns*2] p; // latent probability: strat + outcome
}

transformed parameters {
  // expected counts p[s*y] -> kys[y][s]
  array[2] vector[Ns] kys;
  kys[1] = N * p[0*Ns+1:Ns*1]; // 1st half
  kys[2] = N * p[1*Ns+1:Ns*2]; // 2nd half
  // total counts per (exposure,outcome)
  array[Ne,4] real k4;
  for (i in 1:Ne){
    k4[i,1] = sum(kys[1,se0[i]]); // e=0,y=0
    k4[i,2] = sum(kys[1,se1[i]]); // e=1,y=0
    k4[i,3] = sum(kys[2,se0[i]]); // e=0,y=1
    k4[i,4] = sum(kys[2,se1[i]]); // e=1,y=1
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
  vector[Ns] ks  = kys[1] + kys[2];        // counts per strata
  vector[Ns] pys = kys[2] ./ ks;           // outcome prev per strata
  real pym = sum(ks .* pys) / N;           // weighted mean prevalence
  real pyv = sum(ks .* (pys-pym)^2) / N;   // weighted var prevalence
  real pyd = sqrt(pyv);                    // weighted SD prevalence
  real cv = pyd / pym;                     // coefficient of variation
  array[Ns] int o = sort_indices_asc(pys);     // indices of increasing prev
  vector[Ns] ocps = cumulative_sum(ks[o]) / N; // ord cum prop per strata
  vector[Ns] opys = pys[o];                    // ord prev per strata
}
