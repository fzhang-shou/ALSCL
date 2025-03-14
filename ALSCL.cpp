// age-based length-structured statistical catch-at-length model (ALSCL)
// Z-at-length follows two-way AR1 variation
// decreasing logistic function of the growth increment
// exponential decrease of initial number at age
// configured for quartly time series data

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type pnorm_discrete(Type upper, Type lower){
	
	Type zero = 0.0;
	Type one = 1.0;
	
	Type mid = 0.5 * (upper + lower);
	Type mid2 = mid*mid;
	Type mid4 = mid*mid*mid*mid;
	Type ub = upper - mid;
	Type lb = lower - mid;
	Type d1 = upper - lower;
	Type d3 = ub*ub*ub - lb*lb*lb;
	Type d5 = ub*ub*ub*ub*ub - lb*lb*lb*lb*lb;
	Type ts = d1 + (mid2-1)*d3/6 + (mid4-6*mid2+3)*d5/120; // taylor series expansion
	
	Type ret = dnorm(mid,zero,one,false) * ts;
	
	return ret;
}

template<class Type>
Type objective_function<Type>::operator()()
{
	// data section
	DATA_MATRIX(logN_at_len); // survey number at length
	DATA_MATRIX(na_matrix); // indicating elements where N_at_len is NA.
	DATA_VECTOR(log_q); // survey catchability at length
	DATA_VECTOR(len_mid); // middle value of length bins
	DATA_VECTOR(len_lower); // the lower bound of length bins
	DATA_VECTOR(len_upper); // the upper bound of length bins
	DATA_VECTOR(len_border); // the intervals among length bins (excluding the lower and upper limit)
	DATA_VECTOR(age); // ages of the stock
	DATA_INTEGER(Y);
	DATA_INTEGER(A);
	DATA_INTEGER(L);
	DATA_MATRIX(weight); // weight at length (L,Y)
	DATA_MATRIX(mat); // maturity at length (L,Y)
	DATA_SCALAR(M); // assumed constant natural mortality for large fish.
	DATA_SCALAR(growth_step); // the growth step of the population (usually annual, sometimes quartly)
	
	// parameter section
	PARAMETER(log_init_Z); // intital total mortality used to set initial age structure
	PARAMETER(log_sigma_log_N0); // variability of initial number at age
	
	PARAMETER(mean_log_R);
	PARAMETER(log_sigma_log_R);
	PARAMETER(logit_log_R);
	
	PARAMETER(mean_log_F);
	PARAMETER(log_sigma_log_F);
	PARAMETER(logit_log_F_y);
	PARAMETER(logit_log_F_l);
	
	PARAMETER(log_vbk);
	PARAMETER(log_Linf);
	PARAMETER(log_t0);
	PARAMETER(log_cv_len); // coefficient of variance of the length distribution at age in the first year
	PARAMETER(log_cv_grow); // the coefficient of vairance of growth increment
	PARAMETER(log_sigma_index);
	
	PARAMETER_VECTOR(dev_log_R); // year effect of R as random effects
	PARAMETER_ARRAY(dev_log_F); // deviation of Z across lengths and years
	PARAMETER_VECTOR(dev_log_N0); // year effect of N0 as random effects
	
	// derived parameters
	Type one = 1.0;
	Type zero = 0.0;
	
	Type init_Z = exp(log_init_Z);
	Type sigma_log_N0 = exp(log_sigma_log_N0);
	
	Type sigma_log_R = exp(log_sigma_log_R);
	Type phi_log_R = exp(logit_log_R)/(one+exp(logit_log_R));
	
	Type sigma_log_F = exp(log_sigma_log_F);
	Type phi_log_F_y = exp(logit_log_F_y)/(one+exp(logit_log_F_y));
	Type phi_log_F_l = exp(logit_log_F_l)/(one+exp(logit_log_F_l));
	
	Type vbk = exp(log_vbk);
	Type Linf = exp(log_Linf);
	Type t0 = exp(log_t0);
	Type cv_len = exp(log_cv_len);
	Type cv_grow = exp(log_cv_grow);
	Type sigma_index = exp(log_sigma_index);
	
	using namespace density;
	
	matrix<Type> ZL(L,Y);
	matrix<Type> FL(L,Y);
	matrix<Type> log_FL(L,Y);
	for(int i=0; i<L; ++i){
		for(int j=0; j<Y; ++j){
			log_FL(i,j) = mean_log_F + dev_log_F(i,j);
			FL(i,j) = exp(log_FL(i,j));
			ZL(i,j) = FL(i,j) + M;
		}
	}
	
	// growth transition matrix
	matrix<Type> G(L,L);
	for(int j=0; j<L; ++j){
		
		// Type mean_inc = (one - exp(-vbk))*(Linf - len_mid(j));
		
		Type max_inc = (one - exp(-vbk)) * Linf;
		Type l50 = 0.5 * Linf;
		Type l95 = 0.05 * Linf;
		Type mean_inc = growth_step * max_inc/(one + exp(-log(19.0)*(len_mid(j)-l50)/(l95-l50)));
		
		Type sd_inc = mean_inc * cv_grow;
				
		for(int i=0; i<L; ++i){
			if(i<j){
				G(i,j) = zero; // shrinkage of body length is not allowed
			}
			if(i==j){
				Type mu_lower = len_upper(i)-len_mid(j);
				Type std_mu_lower = (mu_lower-mean_inc)/sd_inc;
				
				int expand = 50; // the expansion of lower and upper bounds to increase the accuracy of estimation for lower and upper bins
				vector<Type> lower_bin(expand);
				Type Ub = std_mu_lower;
				Type Lb = Ub - 1.0;
				lower_bin(0) = pnorm_discrete(Ub,Lb);
				for(int n=1; n<expand; ++n){
					Ub = Ub - 1.0;
					Lb = Ub - 1.0;
					lower_bin(n) = pnorm_discrete(Ub,Lb);
				}
				G(i,j) = lower_bin.sum(); // probability of staying in current length bin
			}
			if(i>j){
				Type mu1 = len_upper(i)-len_mid(j);
			    Type mu2 = len_lower(i)-len_mid(j);
				Type std_mu1 = (mu1-mean_inc)/sd_inc;
				Type std_mu2 = (mu2-mean_inc)/sd_inc;
				//G(i,j) = pnorm(std_mu1) - pnorm(std_mu2);
				
				int expand = 50; // the expansion of lower and upper bounds to increase the accuracy of estimation for lower and upper bins
				vector<Type> mu_range(expand);
				vector<Type> G_range(expand-1);
				mu_range(0) = std_mu2;
				Type inc_mu = (std_mu1 - std_mu2)/(expand-1);
				for(int n=1; n<expand; ++n){
					mu_range(n) = mu_range(n-1) + inc_mu;
					G_range(n-1) = pnorm_discrete(mu_range(n),mu_range(n-1));
				}
				G(i,j)= G_range.sum(); // probability of increasing from length bin j to i
			}
		}
	}
	G((L-1),(L-1)) = one; // probability of staying in the largest length bin is one
	
	for(int i=0; i<L; ++i){
		for(int j=0; j<L; ++j){
			if(G(i,j)<1e-20){G(i,j)=1e-20;};
		}
	}
	
	// get annual recruitment
	vector<Type> log_Rec = dev_log_R + mean_log_R; // may need to plug in SR model here
	vector<Type> Rec = exp(log_Rec); // - half * std_log_R * std_log_R); // bias correction
 
	// length distribution at age for recruitemnt and initial year
	matrix<Type> pla(L,A); // probability of length at age
	for(int j = 0;j < A;++j){
		Type ml = Linf*(one - exp(-vbk*(age(j)-t0)));
		Type sl = ml*cv_len;
		vector<Type> len_border_std = (len_border-ml)/sl; // standardized border among length bins
		
		int expand = 10; // the expansion of lower and upper bins
		vector<Type> lower_bin(expand);
		Type Ub = len_border_std(0);
		Type Lb = Ub - 1.0;
		lower_bin(0) = pnorm_discrete(Ub,Lb);
		for(int n=1; n<expand; ++n){
			Ub = Ub - 1.0;
			Lb = Ub - 1.0;
			lower_bin(n) = pnorm_discrete(Ub,Lb);
		}
		pla(0,j) = lower_bin.sum(); // first length bin
		
		for(int i=1; i< (L-1);++i){
			pla(i,j) = pnorm_discrete(len_border_std(i),len_border_std(i-1));
		}
		
		vector<Type> upper_bin(expand);
		Lb = len_border_std(L-2);
		Ub = Lb + 1.0;
		upper_bin(0) = pnorm_discrete(Ub,Lb);
		for(int n=1; n<expand; ++n){
			Lb = Lb + 1.0;
			Ub = Lb + 1.0;
			upper_bin(n) = pnorm_discrete(Ub,Lb);
		}
		pla(L-1,j) = upper_bin.sum(); // last length bin		
	}
	
	for(int i=0; i<L; ++i){
		for(int j=0; j<A; ++j){
			if(pla(i,j)<1e-20){pla(i,j)=1e-20;};
		}
	}
	// The cohort model;
	array<Type> NLA(L,A,Y);
	array<Type> log_NLA(L,A,Y);
	
	//initializing recruitment and first year abundance
	vector<Type> log_N0(A);
	vector<Type> N0(A);
	log_N0(0) = log_Rec(0);
	N0(0) = Rec(0);
	
	for(int j = 1;j < A;++j){ 
		log_N0(j) = log_N0(j-1) - init_Z + dev_log_N0(j-1);;
		N0(j) = exp(log_N0(j));
	}
	
	for(int j=0; j<A; ++j){
		NLA.col(0).col(j) = pla.col(j) * N0(j); // number at length for each age in the first year
		log_NLA.col(0).col(j) = log(NLA.col(0).col(j));
	}	
	
	// core dynamics
	for(int i = 1; i < Y;++i){
		NLA.col(i).col(0) = pla.col(0) * Rec(i); // stationary length distribution for recruitment
		log_NLA.col(i).col(0) = log(NLA.col(i).col(0));
		for(int j = 1;j < (A-1);++j){
			vector<Type> previous = log_NLA.col(i-1).col(j-1);
			vector<Type> mortality = ZL.col(i-1);
			vector<Type> log_survival = previous - mortality;
			vector<Type> survival = exp(log_survival);
			NLA.col(i).col(j) = G * survival;
			log_NLA.col(i).col(j) = log(NLA.col(i).col(j));
		}
		// plus group
		vector<Type> previous_1 = log_NLA.col(i-1).col(A-2);
		vector<Type> previous_2 = log_NLA.col(i-1).col(A-1);
		vector<Type> mortality = ZL.col(i-1);
		vector<Type> log_survival_1 = previous_1 - mortality;
		vector<Type> log_survival_2 = previous_2 - mortality;
		vector<Type> survival = exp(log_survival_1) + exp(log_survival_2);
		NLA.col(i).col(A-1) = G * survival;
		log_NLA.col(i).col(A-1) = log(NLA.col(i).col(A-1));
	}
	
	// derive three-dimensional total and spawning biomass at length and age
	array<Type> BLA(L,A,Y);
	array<Type> SBLA(L,A,Y);
	for(int i=0; i<Y; ++i){
		for(int j=0; j<A; ++j){
			for(int k=0; k<L; ++k){
				BLA(k,j,i) = NLA(k,j,i)*weight(k,i);
				SBLA(k,j,i) = BLA(k,j,i)*mat(k,i); 				
			}
		}
	}
	
	// derive number, biomass and spawning biomass at length and number at age
	matrix<Type> NL(L,Y); 
	matrix<Type> NA(A,Y);
	matrix<Type> BL(L,Y);
	matrix<Type> BA(A,Y);
	matrix<Type> SBL(L,Y);
	matrix<Type> SBA(A,Y);
	
	for(int i=0; i<Y; ++i){
		for(int j=0; j<A; ++j){
			for(int k=0; k<L; ++k){
				NL(k,i) += NLA(k,j,i); 
				BL(k,i) +=	BLA(k,j,i);
				SBL(k,i) +=SBLA(k,j,i); 				
			}
		}
	}
	for(int i=0; i<Y; ++i){
		for(int j=0; j<A; ++j){
			NA(j,i) = NLA.col(i).col(j).sum();
			BA(j,i) = BLA.col(i).col(j).sum();
			SBA(j,i) = SBLA.col(i).col(j).sum();
		}
	}
	
	// derive catch statistics
	array<Type> CNLA(L,A,Y);
	matrix<Type> CNL(L,Y);
	matrix<Type> CNA(A,Y);
	
	array<Type> CBLA(L,A,Y);
	matrix<Type> CBL(L,Y);
	matrix<Type> CBA(A,Y);
	
	for(int i=0; i<Y; ++i){
		for(int j=0; j<A; ++j){
			for(int k=0; k<L; ++k){
				CNLA(k,j,i) = NLA(k,j,i) * (one - exp(-ZL(k,i))) * (FL(k,i)/(ZL(k,i)));
				CBLA(k,j,i) = CNLA(k,j,i)*weight(k,i);
			}
		}
	}
	
	for(int i=0; i<Y; ++i){
		for(int j=0; j<A; ++j){
			for(int k=0; k<L; ++k){
				CNL(k,i) += CNLA(k,j,i); 
				CBL(k,i) += CBLA(k,j,i); 
			}
		}
	}
	
	for(int i=0; i<Y; ++i){
		for(int j=0; j<A; ++j){
			CNA(j,i) = CNLA.col(i).col(j).sum();
			CBA(j,i) = CBLA.col(i).col(j).sum();
		}
	}
	
	// F, M and Z at age as emergent property derived from population abundance at age and catch at age
	matrix<Type> ZA(A-1,Y-1);
	matrix<Type> FA(A-1,Y-1);
	for(int i=0; i<A-1; ++i){
		for(int j=0; j<Y-1; ++j){
			ZA(i,j) = log(NA(i,j))-log(NA(i+1,j+1));
			FA(i,j) = ZA(i,j) - M;
		}
	}
	
	// derive total abundance, biomass, spawning stock biomass and catch
	vector<Type> N(Y);
	vector<Type> B(Y);
	vector<Type> SSB(Y);
	vector<Type> CN(Y);
	vector<Type> CB(Y);
	for(int i=0; i<Y; ++i){
		N(i) = NL.col(i).sum();
		B(i) = BL.col(i).sum();
		SSB(i) = SBL.col(i).sum();
		CN(i) = CNL.col(i).sum();
		CB(i) = CBL.col(i).sum();
	}
	
	// derive survey index
	matrix<Type> Elog_index(L,Y);
	matrix<Type> resid_index(L,Y);
	for(int i = 0;i < L;++i){
		for(int j=0; j<Y; ++j){
		Elog_index(i,j) = log_q(i) + log(NL(i,j));
		resid_index(i,j) = (logN_at_len(i,j) - Elog_index(i,j))* na_matrix(i,j);
		}
	}
	
	// objective function
	Type nll = zero;
	
	//index, the measurement error
	for(int i=0;i<L;++i){
		for(int j=0;j<Y;++j){
                	if(na_matrix(i,j)>0){
				nll -= dnorm(resid_index(i,j),zero,sigma_index,true);
			}
		}
	}
  
	// logF
	nll += SCALE(SEPARABLE(AR1(phi_log_F_y),AR1(phi_log_F_l)),sigma_log_F)(dev_log_F); // two-way AR1 variation

	//recruitment estimates
	nll += SCALE(AR1(phi_log_R),sigma_log_R)(dev_log_R); // AR1 evaluation to generate nll for logR
	
	// initial population size
	nll -= dnorm(dev_log_N0,zero,sigma_log_N0,true).sum(); // iid evaluation of the initial population size
	
	// report results
	REPORT(pla);
	REPORT(G);
	
	REPORT(Rec);
	REPORT(N);
	REPORT(NLA)
	REPORT(NL);
	REPORT(NA);
	
	REPORT(B);
	REPORT(BLA)
	REPORT(BL);
	REPORT(BA);
	
	REPORT(SSB);
	REPORT(SBLA)
	REPORT(SBL);
	REPORT(SBA);
	
	REPORT(CN)
	REPORT(CNLA);
	REPORT(CNL);
	REPORT(CNA);
	REPORT(CB)
	REPORT(CBLA);
	REPORT(CBL);
	REPORT(CBA);
	
	REPORT(FL);
	REPORT(FA);
	REPORT(ZL);
	REPORT(ZA);
	
	REPORT(mean_log_F);
	REPORT(dev_log_F);
	REPORT(phi_log_F_y);
	REPORT(phi_log_F_l);
	REPORT(sigma_log_F);
	
	REPORT(mean_log_R);
	REPORT(dev_log_R);
	REPORT(sigma_log_R);
	REPORT(phi_log_R);
	
	REPORT(sigma_index);
	REPORT(Elog_index); 
	REPORT(resid_index);
	
	REPORT(vbk);
	REPORT(Linf);
	REPORT(t0);
	REPORT(cv_len);
	REPORT(cv_grow);
	
	REPORT(init_Z);
	REPORT(sigma_log_N0);

    ADREPORT(Rec);
    ADREPORT(N);
    ADREPORT(NA);
    ADREPORT(NL);
    ADREPORT(B);
    ADREPORT(BA);
    ADREPORT(BL);
    ADREPORT(SSB);
    ADREPORT(SBA);
    ADREPORT(SBL);
    ADREPORT(FA);
    ADREPORT(FL);
    ADREPORT(CN);
    ADREPORT(CNA);
    ADREPORT(CNL);

	
	return nll;
}

