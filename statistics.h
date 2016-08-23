/*
 *  RTI (real-time inference)
 *  Statistics.cpp, the header file for all the statistical methods
 *  Currently, maximum likelihood ratio, kolmogorov-smirnov, chi-squared, G, and an AIC method are implemented.
 *  More statistics will be added later (for instance an L-infinity norm bound like many PAC learners).
 *
 *  The sizes of the TIME histogram bins used to calculate the statistics are set to the Interquartile ranges
 *  These should work fine for most time distributions.
 *  The functions in this file are all called from timed_automaton.cpp (merge and split tests)
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

#ifndef _STATISTICS_H_
#define _STATISTICS_H_

class state_statistics;

#include "timed_data.h"
#include "timed_automaton.h"

extern double MAX_DIST;
extern int MIN_DATA;
extern double MAX_P_VALUE;
extern double MIN_P_VALUE;

extern void initialize_consensus_test();
extern void add_to_consensus_test(double p_value);
extern double calculate_consensus_test();
	
extern void initialize_likelihood_test();
extern void add_to_likelihood_test(double ratio, int parameters);
extern double calculate_likelihood_test();
	
extern double calculate_chi2_score(timed_state* old_target, timed_state* new_target);
extern double calculate_chi2_score(timed_state* target);
extern double calculate_chi2_score_time(timed_state* old_target, timed_state* new_target);
extern double calculate_chi2_score_time(timed_state* target);
extern pair<int, double> get_likelihood_ratio(timed_state* old_target, timed_state* new_target);
extern pair<int, double> get_likelihood_ratio(timed_state* target);
extern pair<int, double> get_likelihood_ratio_time(timed_state* old_target, timed_state* new_target);
extern pair<int, double> get_likelihood_ratio_time(timed_state* target);

class state_statistics{
	int total_counts;
	int* symbol_counts;
	int* time_counts;
	
	int total_marks;
	int* symbol_marks;
	int* time_marks;

	friend void initialize_consensus_test();
	friend void add_to_consensus_test(double p_value);
	friend double calculate_consensus_test();
	
	friend void initialize_likelihood_test();
	friend void add_to_likelihood_test(double ratio, int parameters);
	friend double calculate_likelihood_test();
	
	friend double calculate_chi2_score(timed_state* old_target, timed_state* new_target);
	friend double calculate_chi2_score(timed_state* target);
	friend double calculate_chi2_score_time(timed_state* old_target, timed_state* new_target);
	friend double calculate_chi2_score_time(timed_state* target);
	friend pair<int, double> get_likelihood_ratio(timed_state* old_target, timed_state* new_target);
	friend pair<int, double> get_likelihood_ratio(timed_state* target);
	friend pair<int, double> get_likelihood_ratio_time(timed_state* old_target, timed_state* new_target);
	friend pair<int, double> get_likelihood_ratio_time(timed_state* target);

public:	
	state_statistics();
	~state_statistics();
	
	inline int get_time_counts(int t){
		return time_counts[t];
	};
	
	inline int get_symbol_counts(int s){
		return symbol_counts[s];
	};

	const inline int get_bar(int time){
		if(time <= TIME_IQR25) return 0;
		if(time <= TIME_IQR50) return 1;
		if(time <= TIME_IQR75) return 2;
		return 3;
	};
	
	const inline int get_begin_time(int bar){
		switch(bar){
			case 0 : return 0;
			case 1 : return TIME_IQR25 + 1;
			case 2 : return TIME_IQR50 + 1;
			default : return TIME_IQR75 + 1;
		}
	};
	
	const inline int get_end_time(int bar){
		switch(bar){
			case 0 : return TIME_IQR25;
			case 1 : return TIME_IQR50;
			case 2 : return TIME_IQR75;
			default : return MAX_TIME + 1;
		}
	};

	void add_count(timed_tail* tail);
	void del_count(timed_tail* tail);
	void mark(timed_tail* tail);
	void unmark(timed_tail* tail);
	
	inline void add_count(int symbol, int time){
		symbol_counts[symbol]++;
		time_counts[time]++;
		total_counts++;
	};
	
	inline double get_probability(int symbol, int time){
		return (((double) symbol_counts[symbol] * time_counts[time]) / ((double)total_counts * total_counts));
	};

	inline double get_probability_time(int symbol, int time){
		double count = ((double)total_counts / 1000.0) + (double)symbol_counts[symbol];
		double timec = ((double)total_counts / 1000.0) + (double)time_counts[get_bar(time)];
		double additional_count = ((double)total_counts / 1000.0) * (double)MAX_SYMBOL;
		double additional_time = ((double)total_counts / 1000.0) * (double)NUM_HISTOGRAM_BARS;
		return ((count * timec) / (((double)total_counts + additional_count) * ((double)total_counts + additional_time)));
	};

	inline void clear_marks(){
		total_marks = 0;
		for(int i = 0; i < MAX_SYMBOL; ++i){
			symbol_counts[i] += symbol_marks[i];
			symbol_marks[i] = 0;
		}
		
		for(int j = 0; j < NUM_HISTOGRAM_BARS; ++j){
			time_counts[j] += time_marks[j];
			time_marks[j] = 0;
		}
	};

	inline int get_total_counts(){
		return total_counts;
	};

	inline int get_total_marks(){
		return total_marks;
	};

	double get_probability(timed_tail* tail);	
	double get_mark_probability(timed_tail* tail);
};

#endif /* _STATISTICS_H_ */
