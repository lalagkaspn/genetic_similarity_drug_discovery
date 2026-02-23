data{
  //https://mc-stan.org/docs/2_21/stan-users-guide/ragged-data-structs-section.html
  ///
	int <lower=1> N_di; // tot
	int <lower=1> N_top; // tot
	int <lower=1> N_other; // tot
	int <lower=1> N_label; // tot
	int <lower=1> N_evidence;
	int evidence_to_di_top[N_top]; // the non-top
	int evidence_to_di_other[N_other]; // the non-top
	
  int num_per_evidences_top[N_evidence];
  int num_per_evidences_other[N_evidence];
	int<lower=0, upper=1> indicated[N_label];
	int<lower=1, upper=N_di> drug_ind_index_label[N_label]; // matching "indicated" labels to the drug jndications
	vector[N_top] cors_top;
	vector[N_other] cors_other;
	//vector[N] body_system;
	//real cor_sd_wid_top;
	real sd_multiplier;
	vector[N_evidence] cor_sd_manual;
	//real cor_sd_wid_other;
  real<lower=0,upper=1> phi;
  real<lower=0.1> lambda;
  //vector<lower=0,upper=1> phiTheta[N_evidence];
  //vector<lower=0.1> lambdaTheta[N_evidence];
  //real cor_sd_mean;
  int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood
  //int<lower = 0, upper = 1> hurdle_zeros[N_evidence];
}
parameters{
	//real<lower=0, upper=1> p_di[N_di];
	vector<lower=0, upper=1>[N_di] p_di;
	vector<lower=0, upper=1>[N_evidence] intercept_cor; //
	row_vector<lower=0, upper=1>[N_evidence] slope_cor; // note - removing hte bounds on this leads to divergences??
	//real<lower=0> cor_sd_top;
	//real<lower=0> cor_sd_other;
	//real<lower=0, upper=1> intercept_indicated; // change this for intercept per body so 12? intercepts
  //real<lower=0, upper=1> theta;
   // real<lower=1, upper=2> sd_multiplier;
  
}
transformed parameters{
  //real<lower=0, upper=1> prob_for_sampling[N_di]; // = p_di + rep_array(intercept_indicated, N_di); // intercept_indicated; 
  //https://mc-stan.org/docs/2_18/reference-manual/vector-and-matrix-data-types.html
  //vector<lower=0, upper=1>[N] prob_for_sampling;
   vector<lower=0, upper=1>[N_di] cor_means[N_evidence];//
  //vector<lower=0, upper=1>[N_di] cor_for_sampling_top;//
  //real log_theta = log(theta);
  //real log_1mtheta = log1m(theta);

	
	//for (di in 1:N_di){
	//  prob_for_sampling[di] = p_di[di]; // + intercept_indicated;
	//}
	for(n in 1:N_evidence){
	  cor_means[n] = intercept_cor[n] + slope_cor[n] * p_di;
	}
	//cor_means = slope_cor * rep_matrix(p_di, N_evidence); //intercept_cor//intercept_cor_top[body_system] //[drug_ind_index]; 
  //cor_for_sampling_other = intercept_cor + slope_cor * p_di; //[drug_ind_index]; 

	//for (n in 1:N){	
	     //prob_for_sampling[n] = p_di[drug_ind_index[n]] ;//+ intercept_indicated;
	 //    cor_for_sampling[n] = p_di[drug_ind_index[n]]*slope_cor + intercept_cor;
	 //}
  //cor_for_sampling = intercept_cor  + slope_cor * prob_for_sampling; //, cor_var
}

model {
    //https://en.wikipedia.org/wiki/Beta_distribution --> "Mean and sample size" parameterization
    real alpha = lambda * phi;
    real beta = lambda * (1 - phi);
    //real alphaTheta = lambdaTheta * phiTheta;
    //real betaTheta = lambdaTheta * (1 - phiTheta);
   // cor_sd_top ~ normal(cor_sd_mean, cor_sd_wid_top);
    //cor_sd_other ~ normal(0, cor_sd_wid_other);
    p_di ~ beta(alpha, beta);
    //theta ~ beta(alphaTheta, betaTheta);
    slope_cor ~ normal(0, .2);
    intercept_cor ~ normal(0,.2);
    
    indicated ~ bernoulli(p_di[drug_ind_index_label]);
    int pos_top = 1;
    int pos_other = 1;
    // segment(cors_top, 1, num_per_evidences_top[1]) ~normal(cor_means[1][segment(evidence_to_di_top, 1, num_per_evidences_top[1])],
    //                                                             cor_sd_manual[1]); // COMPILES
    for(evidence in 1:N_evidence){
      segment(cors_top, pos_top, num_per_evidences_top[evidence]) ~ normal(cor_means[evidence][segment(evidence_to_di_top, pos_top, num_per_evidences_top[evidence])],
                                                                cor_sd_manual[evidence]);
      pos_top = pos_top + num_per_evidences_top[evidence];
      segment(cors_other, pos_other, num_per_evidences_other[evidence]) ~ normal(cor_means[evidence][segment(evidence_to_di_other, pos_other, num_per_evidences_other[evidence])],
                                                                cor_sd_manual[evidence]*sd_multiplier);
      pos_other = pos_other + num_per_evidences_other[evidence];
    }
}
      // 	      for (n in 1:N){
//           if (cors[n] != 0) {
//              real lpdf = normal_lpdf(cors[n] | cor_for_sampling[drug_ind_index[n]], cor_sd);
//              target += lpdf;
//           }else{
//            n_zero += 1 ;
//           }
// 	      }
//         target += binomial_lpmf(n_zero | N, theta);

	 //     target += normal_lpdf(cors_top | cor_for_sampling_top, cor_sd_manual);
  	 //     target += normal_lpdf(cors_other | cor_for_sampling_other[drug_ind_index], sd_multiplier*cor_sd_manual);
    
//     if(run_estimation==1){
//       indicated ~ bernoulli_logit(p_di); //prob_for_sampling); //intercept_indicated + prob_for_sampling);
// 	    if(hurdle_zeros==1){
// 	      int n_zero = 0;
// 	      for (n in 1:N){
//           if (cors[n] != 0) {
//              real lpdf = normal_lpdf(cors[n] | cor_for_sampling[drug_ind_index[n]], cor_sd);
//              target += lpdf;
//           }else{
//            n_zero += 1 ;
//           }
// 	      }
//         target += binomial_lpmf(n_zero | N, theta);
// 	    }else{
// 	      target += normal_lpdf(cors_top, cor_for_sampling_top, cor_sd_top);
// 	      target += normal_lpdf(cors_other, cor_for_sampling[drug_ind_index], cor_sd_other);
// 	    }
// 
// 	 }
 
 
// generated quantities{
// 	real<lower=0, upper=1> p_di_sim[N_di];
// 
// 	real<lower=0, upper=1> cor_for_sampling_sim[N];
// 	  int<lower = 0, upper = 1> indicated_sim[N_di];
// 
//  	  vector[N] cor_rep2; //<lower=0, upper=1>
// 	  vector[2] theta_samp; //I think this needs to have
// 	  int<lower = 1, upper = 2> notzero[N];
// 	  theta_samp[1] = 1 - theta;
// 	  theta_samp[2] = theta;
// 	 p_di_sim = beta_rng(rep_array(lambda * phi, N_di),
// 	                    rep_array(lambda * (1 - phi),N_di));
//    indicated_sim = bernoulli_rng(p_di_sim);
//    if(hurdle_zeros==0){
//      for (n in 1:N) {
//        notzero[n]= 1;
//       cor_for_sampling_sim[n] = intercept_cor + slope_cor * p_di_sim[drug_ind_index[n]]; //rep_array(intercept_cor,N) + rep_array(slope_cor,N) * p_di_sim[drug_ind_index] ; //* rep_cor; //rep_array(slope_cor,N) .*
//       cor_rep2[n] = normal_rng(cor_for_sampling_sim[n], cor_sd);
//      }
//    }else{
//      for (n in 1:N) {
//        cor_for_sampling_sim[n] = intercept_cor + slope_cor * p_di_sim[drug_ind_index[n]];
//        notzero[n]  = categorical_rng(theta_samp);
//        if(notzero[n] == 1){
//          cor_rep2[n] = 0;
//        }else{ //==2
//          cor_rep2[n] = normal_rng(cor_for_sampling_sim[n], cor_sd);
//        }
//     }
//    }
// 
// }


