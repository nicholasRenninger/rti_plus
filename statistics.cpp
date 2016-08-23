/*
 *  RTI (real-time inference)
 *  Statistics.cpp, the source file for all the statistical methods
 *  Currently, maximum likelihood ratio, and the chi-squared Fisher's method (best performing) from my PhD thesis are implemented.
 *  More statistics will be added later (for instance an L-infinity norm bound like many PAC learners).
 *
 *  The functions ar all called from timed_automaton.cpp
 *
 *  For some info regarding the algorithm, see:
 *  Sicco Verwer and Mathijs de Weerdt and Cees Witteveen (2007),
 *  An algorithm for learning real-time automata,
 *  In Maarten van Someren and Sophia Katrenko and Pieter Adriaans (Eds.),
 *  Proceedings of the Sixteenth Annual Machine Learning Conference of Belgium and the Netherlands (Benelearn),
 *  pp. 128-135.
 *  
 *  Copyright 2009 - Sicco Verwer, jan-2009
 *  This program is released under the GNU General Public License
 *  Info online: http://www.gnu.org/licenses/quick-guide-gplv3.html
 *  Or in the file: licence.txt
 *  For information/questions contact: siccoverwer@gmail.com
 *
 *  I will try to keep this software updated and will also try to add new statistics or search methods.
 *  Also I will add comments to the source in the near future.
 *
 *  Feel free to adapt the code to your needs, please inform me of (potential) improvements.
 */

#include <math.h>
#include <gsl/gsl_cdf.h>
#include "statistics.h"

double MAX_DIST = 0.05;
int MIN_DATA = 10;
double MAX_P_VALUE = 1.0 - 0.1e-100;
double MIN_P_VALUE = 0.1e-100;

/* Calculates the G test */
double calculate_G_value(double first, double second, double total1, double total2){
	double total     = (double)(first + second);
	double expected1 = (total1 * total) / (total1 + total2);
	double expected2 = (total2 * total) / (total1 + total2);
	
	return 2.0 * first * log(first / expected1) + 2.0 * second * log(second/expected2);
};

/* Calculates the chi^2 test */
double calculate_chi2_value(double first, double second, double total1, double total2){
	double total     = (double)(first + second);
	double expected1 = (total1 * total) / (total1 + total2);
	double expected2 = (total2 * total) / (total1 + total2);

	double top1 = first - expected1;
	double top2 = second - expected2;
	
	/* Yates correction for continuity */
	if(first < MIN_DATA || second < MIN_DATA){
		if(top1 < 0) top1 = -top1;
		top1 -= 0.5;
		if(top2 < 0) top2 = -top2;
		top2 -= 0.5;
	}
	
	return ((top1*top1)/expected1) + ((top2*top2)/expected2);
};

/* The Fisher's method consensus test */
double sum_z_values = 0.0;
double num_tests = 0.0;

void initialize_consensus_test(){
	sum_z_values = 0.0;
	num_tests = 0.0;
};

void add_to_consensus_test(double p_value){
	if(p_value == 1.0) p_value = MAX_P_VALUE;
	sum_z_values += - 2.0 * log(p_value);
	num_tests++;
};

double calculate_consensus_test(){
	if(num_tests == 0) return -1.0;
	return gsl_cdf_chisq_Q (sum_z_values, (double)(2 * num_tests));
};
/* End of Fisher's method consensus test */

/* The likelihood ratio test */
double ml_ratio = 0.0;
int ml_parameters = 0;

void initialize_likelihood_test(){
	ml_ratio = 0.0;
	ml_parameters = 0;
};

void add_to_likelihood_test(double ratio, int parameters){
	ml_ratio += ratio;
	ml_parameters += parameters;
};

double calculate_likelihood_test(){
	if(ml_ratio == 0.0 && ml_parameters == 0) return -1.0;
	double chi2_value = - 2.0 * ml_ratio;
	double p_value = gsl_cdf_chisq_Q (chi2_value, ml_parameters);
	return p_value;
};
/* End of likelihood ratio test */

/* Calculates the chi^2 value of the SYMBOL distributions
 * for merging two states and adds the result to the consensus test */ 
double calculate_chi2_score(timed_state* old_target, timed_state* new_target){
	if(old_target == 0 || new_target == 0) return -1.0;
	
	/* total counts */
	int total_old = old_target->stat->get_total_counts();
	int total_new = new_target->stat->get_total_counts();
	if(total_old < MIN_DATA || total_new < MIN_DATA) return -1.0;

	/* pooling less than MIN_DATA counts */
	int old_pool = 0;
	int new_pool = 0;
	for(int i = 0; i < MAX_SYMBOL; ++i){
		if(old_target->stat->symbol_counts[i] < MIN_DATA && new_target->stat->symbol_counts[i] < MIN_DATA){
			old_pool += old_target->stat->symbol_counts[i];
			new_pool += new_target->stat->symbol_counts[i];
		}
	}
	
	if(old_pool < MIN_DATA && new_pool < MIN_DATA){
		total_old -= old_pool;
		total_new -= new_pool;
		old_pool = 0;
		new_pool = 0;
	}
	
	/* calculating chi2 dof (degree of freedom) and value */
	double chi2_value = 0.0;
	double chi2_dof = -1.0;
	
	for(int i = 0; i < MAX_SYMBOL; ++i){
		if(old_target->stat->symbol_counts[i] < MIN_DATA && new_target->stat->symbol_counts[i] < MIN_DATA) continue;
			
		chi2_value += calculate_chi2_value(old_target->stat->symbol_counts[i], new_target->stat->symbol_counts[i], total_old, total_new);
		chi2_dof += 1.0;
	}
	
	if(old_pool > MIN_DATA || new_pool > MIN_DATA){
		chi2_value += calculate_chi2_value(old_pool, new_pool, total_old, total_new);
		chi2_dof  += 1.0;
	}
	
	/* testing and adding to consensus test */
	if(chi2_dof >= 1.0){
		double p_value = gsl_cdf_chisq_Q (chi2_value, chi2_dof);
		
		if(p_value < MIN_P_VALUE) p_value = MIN_P_VALUE;
		add_to_consensus_test(p_value);
		return p_value;
	}
	return -1.0;
};

/* Calculates the chi^2 value of the SYMBOL distributions
 * for splitting a state and adds the result to the consensus test */ 
double calculate_chi2_score(timed_state* target){
	if(target == 0) return -1.0;
	
	/* total counts */
	int total_old = target->stat->get_total_counts();
	int total_new = target->stat->get_total_marks();
	if(total_old < MIN_DATA || total_new < MIN_DATA) return -1.0;

	/* pooling less than MIN_DATA counts */
	int old_pool = 0;
	int new_pool = 0;
	for(int i = 0; i < MAX_SYMBOL; ++i){
		if(target->stat->symbol_counts[i] < MIN_DATA && target->stat->symbol_marks[i] < MIN_DATA){
			old_pool += target->stat->symbol_counts[i];
			new_pool += target->stat->symbol_marks[i];
		}
	}
	
	if(old_pool < MIN_DATA && new_pool < MIN_DATA){
		total_old -= old_pool;
		total_new -= new_pool;
		old_pool = 0;
		new_pool = 0;
	};
	
	/* calculating chi2 dof (degree of freedom) and value */
	double chi2_value = 0.0;
	double chi2_dof = -1.0;
	
	for(int i = 0; i < MAX_SYMBOL; ++i){
		if(target->stat->symbol_counts[i] < MIN_DATA && target->stat->symbol_marks[i] < MIN_DATA) continue;
		
		chi2_value += calculate_chi2_value(target->stat->symbol_counts[i], target->stat->symbol_marks[i], total_old, total_new);
		chi2_dof += 1.0;
	}
	
	if(old_pool > MIN_DATA || new_pool > MIN_DATA){
		chi2_value += calculate_chi2_value(old_pool, new_pool, total_old, total_new);
		chi2_dof  += 1.0;
	}
	
	/* testing and adding to consensus test */
	if(chi2_dof >= 1.0){
		double p_value = gsl_cdf_chisq_Q (chi2_value, chi2_dof);
		
		if(p_value < MIN_P_VALUE) p_value = MIN_P_VALUE;
		add_to_consensus_test(p_value);
		return p_value;
	}
	return -1.0;
};

/* Calculates the chi^2 value of the TIME distributions
 * for merging two states and adds the result to the consensus test */ 
double calculate_chi2_score_time(timed_state* old_target, timed_state* new_target){
	if(old_target == 0 || new_target == 0) return -1.0;
	
	/* total counts */
	int total_old = old_target->stat->get_total_counts();
	int total_new = new_target->stat->get_total_counts();
	if(total_old < MIN_DATA || total_new < MIN_DATA) return -1.0;

	/* pooling less than MIN_DATA counts */
	int old_pool = 0;
	int new_pool = 0;
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i){
		if(old_target->stat->time_counts[i] < MIN_DATA && new_target->stat->time_counts[i] < MIN_DATA){
			old_pool += old_target->stat->time_counts[i];
			new_pool += new_target->stat->time_counts[i];
		}
	}
	
	if(old_pool < MIN_DATA && new_pool < MIN_DATA){
		total_old -= old_pool;
		total_new -= new_pool;
		old_pool = 0;
		new_pool = 0;
	};
	
	/* calculating chi2 dof (degree of freedom) and value */
	double chi2_value = 0.0;
	double chi2_dof = -1.0;
	
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i){
		if(old_target->stat->time_counts[i] < MIN_DATA && new_target->stat->time_counts[i] < MIN_DATA) continue;
			
		chi2_value += calculate_chi2_value(old_target->stat->time_counts[i], new_target->stat->time_counts[i], total_old, total_new);
		chi2_dof += 1.0;
	}
	
	if(old_pool > MIN_DATA || new_pool > MIN_DATA){
		chi2_value += calculate_chi2_value(old_pool, new_pool, total_old, total_new);
		chi2_dof  += 1.0;
	}
	
	/* testing and adding to consensus test */
	if(chi2_dof >= 1.0){
		double p_value = gsl_cdf_chisq_Q (chi2_value, chi2_dof);
		
		if(p_value < MIN_P_VALUE) p_value = MIN_P_VALUE;
		add_to_consensus_test(p_value);
		return p_value;
	}
	return -1.0;
};

/* Calculates the chi^2 value of the TIME distributions
 * for splitting a state and adds the result to the consensus test */ 
double calculate_chi2_score_time(timed_state* target){
	if(target == 0) return -1.0;
	
	/* total counts */
	int total_old = target->stat->get_total_counts();
	int total_new = target->stat->get_total_marks();
	if(total_old < MIN_DATA || total_new < MIN_DATA) return -1.0;

	/* pooling less than MIN_DATA counts */
	int old_pool = 0;
	int new_pool = 0;
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i){
		if(target->stat->time_counts[i] < MIN_DATA && target->stat->time_marks[i] < MIN_DATA){
			old_pool += target->stat->time_counts[i];
			new_pool += target->stat->time_marks[i];
		}
	}
	
	if(old_pool < MIN_DATA && new_pool < MIN_DATA){
		total_old -= old_pool;
		total_new -= new_pool;
		old_pool = 0;
		new_pool = 0;
	};
	
	/* calculating chi2 dof (degree of freedom) and value */
	double chi2_value = 0.0;
	double chi2_dof = -1.0;
	
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i){
		if(target->stat->time_counts[i] < MIN_DATA && target->stat->time_marks[i] < MIN_DATA) continue;
		
		chi2_value += calculate_chi2_value(target->stat->time_counts[i], target->stat->time_marks[i], total_old, total_new);
		chi2_dof += 1.0;
	}
	
	if(old_pool > MIN_DATA || new_pool > MIN_DATA){
		chi2_value += calculate_chi2_value(old_pool, new_pool, total_old, total_new);
		chi2_dof  += 1.0;
	}
	
	/* testing and adding to consensus test */
	if(chi2_dof >= 1.0){
		double p_value = gsl_cdf_chisq_Q (chi2_value, chi2_dof);
		
		if(p_value < MIN_P_VALUE) p_value = MIN_P_VALUE;
		add_to_consensus_test(p_value);
		return p_value;
	}
	return -1.0;
};

/* Calculates the likelihood ratio of the SYMBOL distributions
 * for merging two states and adds the result to the likelihood ratio test */ 
pair<int, double> get_likelihood_ratio(timed_state* old_target, timed_state* new_target){
	int extra_parameters = 0;
	double ratio = 0.0;
	
	/* total counts */
	int total_old = old_target->stat->get_total_counts();
	int total_new = new_target->stat->get_total_counts();
	if(total_old < MIN_DATA || total_new < MIN_DATA) return pair<int, double>(0, 0.0);
	
	/* pooling less than MIN_DATA counts */
	int old_pool = 0;
	int new_pool = 0;
	for(int i = 0; i < MAX_SYMBOL; ++i){
		if(old_target->stat->symbol_counts[i] < MIN_DATA && new_target->stat->symbol_counts[i] < MIN_DATA){
			old_pool += old_target->stat->symbol_counts[i];
			new_pool += new_target->stat->symbol_counts[i];
		}
	}
	
	if(old_pool < MIN_DATA || new_pool < MIN_DATA){
		total_old -= old_pool;
		total_new -= new_pool;
		old_pool = 0;
		new_pool = 0;
	}
	
	/* calculating ratio and parameters */
	for(int i = 0; i < MAX_SYMBOL; ++i){
		if(old_target->stat->symbol_counts[i] < MIN_DATA && new_target->stat->symbol_counts[i] < MIN_DATA) continue;
		
		double top_probability = ((double)(old_target->stat->symbol_counts[i] + new_target->stat->symbol_counts[i]))
		/ ((double)(total_old + total_new));
		
		double bottom_probability1 = 1.0;
		if(old_target->stat->symbol_counts[i] != 0) bottom_probability1 = ((double)old_target->stat->symbol_counts[i]) / ((double)total_old);
		
		double bottom_probability2 = 1.0;
		if(new_target->stat->symbol_counts[i] != 0) bottom_probability2 = ((double)new_target->stat->symbol_counts[i]) / ((double)total_new);
		
		ratio += (double)old_target->stat->symbol_counts[i] * log(top_probability);
		ratio -= (double)old_target->stat->symbol_counts[i] * log(bottom_probability1);
		ratio += (double)new_target->stat->symbol_counts[i] * log(top_probability);
		ratio -= (double)new_target->stat->symbol_counts[i] * log(bottom_probability2);
		extra_parameters++;
	}
	
	if(old_pool > MIN_DATA || new_pool > MIN_DATA){
		double top_probability = ((double)(old_pool + new_pool)) / ((double)(total_old + total_new));
		double bottom_probability1 = 1.0;
		if(old_pool != 0) bottom_probability1 = ((double)old_pool) / ((double)total_old);
		double bottom_probability2 = 1.0;
		if(new_pool != 0) bottom_probability2 = ((double)new_pool) / ((double)total_new);
		
		ratio += (double)old_pool * log(top_probability);
		ratio -= (double)old_pool * log(bottom_probability1);
		ratio += (double)new_pool * log(top_probability);
		ratio -= (double)new_pool * log(bottom_probability2);
		extra_parameters++;
	}
	if(extra_parameters > 0){
		add_to_likelihood_test(ratio, extra_parameters);
		return pair<int, double>(extra_parameters, ratio);
	}
	return pair<int, double>(0, 0.0);
};

/* Calculates the likelihood ratio of the TIME distributions
 * for merging two states and adds the result to the likelihood ratio test */ 
pair<int, double> get_likelihood_ratio_time(timed_state* old_target, timed_state* new_target){
	int extra_parameters = 0;
	double ratio = 0.0;
	
	/* total counts */
	int total_old = old_target->stat->get_total_counts();
	int total_new = new_target->stat->get_total_counts();
	if(total_old < MIN_DATA || total_new < MIN_DATA) return pair<int, double>(0, 0.0);

	int old_pool = 0;
	int new_pool = 0;
	/* pooling less than MIN_DATA counts */
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i){
		if(old_target->stat->time_counts[i] < MIN_DATA && new_target->stat->time_counts[i] < MIN_DATA){
			old_pool += old_target->stat->time_counts[i];
			new_pool += new_target->stat->time_counts[i];
		}
	}
	
	if(old_pool < MIN_DATA || new_pool < MIN_DATA){
		total_old -= old_pool;
		total_new -= new_pool;
		old_pool = 0;
		new_pool = 0;
	}
	
	/* calculating ratio and parameters */
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i){
		if(old_target->stat->time_counts[i] < MIN_DATA && new_target->stat->time_counts[i] < MIN_DATA) continue;
		
		double top_probability = ((double)(old_target->stat->time_counts[i] + new_target->stat->time_counts[i]))
		/ ((double)(total_old + total_new));
		
		double bottom_probability1 = 1.0;
		if(old_target->stat->time_counts[i] != 0) bottom_probability1 = ((double)old_target->stat->time_counts[i]) / ((double)total_old);
		
		double bottom_probability2 = 1.0;
		if(new_target->stat->time_counts[i] != 0) bottom_probability2 = ((double)new_target->stat->time_counts[i]) / ((double)total_new);
		
		ratio += (double)old_target->stat->time_counts[i] * log(top_probability);
		ratio -= (double)old_target->stat->time_counts[i] * log(bottom_probability1);
		ratio += (double)new_target->stat->time_counts[i] * log(top_probability);
		ratio -= (double)new_target->stat->time_counts[i] * log(bottom_probability2);
		extra_parameters++;
	}
	
	if(old_pool > MIN_DATA || new_pool > MIN_DATA){
		double top_probability = ((double)(old_pool + new_pool)) / ((double)(total_old + total_new));
		double bottom_probability1 = 1.0;
		if(old_pool != 0) bottom_probability1 = ((double)old_pool) / ((double)total_old);
		double bottom_probability2 = 1.0;
		if(new_pool != 0) bottom_probability2 = ((double)new_pool) / ((double)total_new);
		
		ratio += (double)old_pool * log(top_probability);
		ratio -= (double)old_pool * log(bottom_probability1);
		ratio += (double)new_pool * log(top_probability);
		ratio -= (double)new_pool * log(bottom_probability2);
		extra_parameters++;
	}
	if(extra_parameters > 0){
		add_to_likelihood_test(ratio, extra_parameters);
		return pair<int, double>(extra_parameters, ratio);
	}
	return pair<int, double>(0, 0.0);
};

/* Calculates the likelihood ratio of the SYMBOL distributions
 * for splitting a state and adds the result to the likelihood ratio test */ 
pair<int, double> get_likelihood_ratio(timed_state* target){
	int extra_parameters = 0;
	double ratio = 0.0;
	
	/* total counts */
	int total_old = target->stat->get_total_counts();
	int total_new = target->stat->get_total_marks();
	if(total_old < MIN_DATA || total_new < MIN_DATA) return pair<int, double>(0, 0.0);;
	
	/* pooling less than MIN_DATA counts */
	int old_pool = 0;
	int new_pool = 0;
	for(int i = 0; i < MAX_SYMBOL; ++i){
		if(target->stat->symbol_counts[i] < MIN_DATA && target->stat->symbol_marks[i] < MIN_DATA){
			old_pool += target->stat->symbol_counts[i];
			new_pool += target->stat->symbol_marks[i];
		}
	}
	
	if(old_pool < MIN_DATA || new_pool < MIN_DATA){
		total_old -= old_pool;
		total_new -= new_pool;
		old_pool = 0;
		new_pool = 0;
	}
	
	/* calculating ratio and parameters */
	for(int i = 0; i < MAX_SYMBOL; ++i){
		if(target->stat->symbol_counts[i] < MIN_DATA && target->stat->symbol_marks[i] < MIN_DATA) continue;
		
		double top_probability = ((double)(target->stat->symbol_counts[i] + target->stat->symbol_marks[i]))
		/ ((double)(total_old + total_new));
		
		double bottom_probability1 = 1.0;
		if(target->stat->symbol_counts[i] != 0) bottom_probability1 = ((double)target->stat->symbol_counts[i]) / ((double)total_old);
		
		double bottom_probability2 = 1.0;
		if(target->stat->symbol_marks[i] != 0) bottom_probability2 = ((double)target->stat->symbol_marks[i]) / ((double)total_new);
		
		ratio += (double)target->stat->symbol_counts[i] * log(top_probability);
		ratio -= (double)target->stat->symbol_counts[i] * log(bottom_probability1);
		ratio += (double)target->stat->symbol_marks[i] * log(top_probability);
		ratio -= (double)target->stat->symbol_marks[i] * log(bottom_probability2);
		extra_parameters++;
	}
	
	if(old_pool > MIN_DATA || new_pool > MIN_DATA){
		double top_probability = ((double)(old_pool + new_pool)) / ((double)(total_old + total_new));
		double bottom_probability1 = 1.0;
		if(old_pool != 0) bottom_probability1 = ((double)old_pool) / ((double)total_old);
		double bottom_probability2 = 1.0;
		if(new_pool != 0) bottom_probability2 = ((double)new_pool) / ((double)total_new);
		
		ratio += (double)old_pool * log(top_probability);
		ratio -= (double)old_pool * log(bottom_probability1);
		ratio += (double)new_pool * log(top_probability);
		ratio -= (double)new_pool * log(bottom_probability2);
		extra_parameters++;
	}
	if(extra_parameters > 0){
		add_to_likelihood_test(ratio, extra_parameters);
		return pair<int, double>(extra_parameters, ratio);
	}
	return pair<int, double>(0, 0.0);
};

/* Calculates the likelihood ratio of the TIME distributions
 * for splitting a state and adds the result to the likelihood ratio test */ 
pair<int, double> get_likelihood_ratio_time(timed_state* target){
	int extra_parameters = 0;
	double ratio = 0.0;
	
	/* total counts */
	int total_old = target->stat->get_total_counts();
	int total_new = target->stat->get_total_marks();
	if(total_old < MIN_DATA || total_new < MIN_DATA) return pair<int, double>(0, 0.0);;
	
	/* pooling less than MIN_DATA counts */
	int old_pool = 0;
	int new_pool = 0;
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i){
		if(target->stat->time_counts[i] < MIN_DATA && target->stat->time_marks[i] < MIN_DATA){
			old_pool += target->stat->time_counts[i];
			new_pool += target->stat->time_marks[i];
		}
	}
	
	if(old_pool < MIN_DATA || new_pool < MIN_DATA){
		total_old -= old_pool;
		total_new -= new_pool;
		old_pool = 0;
		new_pool = 0;
	}
	
	/* calculating ratio and parameters */
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i){
		if(target->stat->time_counts[i] < MIN_DATA && target->stat->time_marks[i] < MIN_DATA) continue;
		
		double top_probability = ((double)(target->stat->time_counts[i] + target->stat->time_marks[i]))
		/ ((double)(total_old + total_new));
		
		double bottom_probability1 = 1.0;
		if(target->stat->time_counts[i] != 0) bottom_probability1 = ((double)target->stat->time_counts[i]) / ((double)total_old);
		
		double bottom_probability2 = 1.0;
		if(target->stat->time_marks[i] != 0) bottom_probability2 = ((double)target->stat->time_marks[i]) / ((double)total_new);
		
		ratio += (double)target->stat->time_counts[i] * log(top_probability);
		ratio -= (double)target->stat->time_counts[i] * log(bottom_probability1);
		ratio += (double)target->stat->time_marks[i] * log(top_probability);
		ratio -= (double)target->stat->time_marks[i] * log(bottom_probability2);
		extra_parameters++;
	}
	
	if(old_pool > MIN_DATA || new_pool > MIN_DATA){
		double top_probability = ((double)(old_pool + new_pool)) / ((double)(total_old + total_new));
		double bottom_probability1 = 1.0;
		if(old_pool != 0) bottom_probability1 = ((double)old_pool) / ((double)total_old);
		double bottom_probability2 = 1.0;
		if(new_pool != 0) bottom_probability2 = ((double)new_pool) / ((double)total_new);
		
		ratio += (double)old_pool * log(top_probability);
		ratio -= (double)old_pool * log(bottom_probability1);
		ratio += (double)new_pool * log(top_probability);
		ratio -= (double)new_pool * log(bottom_probability2);
		extra_parameters++;
	}
	if(extra_parameters > 0){
		add_to_likelihood_test(ratio, extra_parameters);
		return pair<int, double>(extra_parameters, ratio);
	}
	return pair<int, double>(0, 0.0);
};

/* Constructor */
state_statistics::state_statistics(){
	total_counts = 0;
	symbol_counts = new int[MAX_SYMBOL];
	for(int i = 0; i < MAX_SYMBOL; ++i) symbol_counts[i] = 0;
	time_counts = new int[NUM_HISTOGRAM_BARS];
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i) time_counts[i] = 0;

	total_marks = 0;
	symbol_marks = new int[MAX_SYMBOL];
	for(int i = 0; i < MAX_SYMBOL; ++i) symbol_marks[i] = 0;
	time_marks = new int[NUM_HISTOGRAM_BARS];
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i) time_marks[i] = 0;
};

/* Destructor */
state_statistics::~state_statistics(){
	delete[] symbol_counts;
	delete[] time_counts;
	delete[] symbol_marks;
	delete[] time_marks;
};

void state_statistics::add_count(timed_tail* tail){
		total_counts++;
		symbol_counts[tail->get_symbol()]++;
		time_counts[get_bar(tail->get_time_value())]++;
};

void state_statistics::del_count(timed_tail* tail){
		total_counts--;
		symbol_counts[tail->get_symbol()]--;
		time_counts[get_bar(tail->get_time_value())]--;
};

void state_statistics::mark(timed_tail* tail){
		int bar_number = get_bar(tail->get_time_value());
		total_marks++;
		symbol_marks[tail->get_symbol()]++;
		time_marks[bar_number]++;
		total_counts--;
		symbol_counts[tail->get_symbol()]--;
		time_counts[bar_number]--;
};

void state_statistics::unmark(timed_tail* tail){
		int bar_number = get_bar(tail->get_time_value());
		total_marks--;
		symbol_marks[tail->get_symbol()]--;
		time_marks[bar_number]--;
		total_counts++;
		symbol_counts[tail->get_symbol()]++;
		time_counts[bar_number]++;
};

double state_statistics::get_probability(timed_tail* tail){
		int bar_number = get_bar(tail->get_time_value());
		return (double)(symbol_counts[tail->get_symbol()] * time_counts[bar_number]) / (double)(total_counts * total_counts);
};

double state_statistics::get_mark_probability(timed_tail* tail){
		int bar_number = get_bar(tail->get_time_value());
		return (double)(symbol_marks[tail->get_symbol()] * time_marks[bar_number]) / (double)(total_marks * total_marks);
};
