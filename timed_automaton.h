/*
 *  RTI (real-time inference)
 *  Timed_automaton.h, the header file for real-time automata, essentially just a set of states.
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

#ifndef TIMED_AUTOMATON_H_
#define TIMED_AUTOMATON_H_

using namespace std;
#include <vector>
#include <sstream>
#include "interval.h"
#include "statistics.h"

extern double MAX_DIST;
extern int MIN_DATA;
extern int TEST_TYPE;

class timed_automaton;
class timed_state;

typedef vector<timed_state*> state_list;
typedef state_list::iterator state_it;

class timed_automaton{
private:
	state_list states;
	timed_state* root;
	timed_input* input;
	
	void check_next_tail(interval* in, timed_tail* tail);
	void recursive_tree_automaton(timed_state*, timed_state*);
	int recursive_total_num_states(timed_state*);

public:
	timed_automaton();
  timed_automaton(timed_input* in);
	~timed_automaton();
	
	void check_consistency();
	
	string to_str();
	string to_str_full();
	void from_file(FILE* file);
	void tree_automaton();
	void garbage_automaton();
	int total_num_states();
	int get_size();

	inline timed_state* get_root(){
		return root;
	};
	
	inline state_list& get_states(){
		return states;
	};
	
	inline void add_state(timed_state* s){
		states.push_back(s);
	};

	inline void del_state(timed_state* s){
		state_list::iterator it = states.end();
		while(it != states.begin()){
			--it;
			if((*it) == s){
				states.erase(it);
				break;
			}
		}
	};
	
	inline bool contains_state(timed_state* s){
		state_list::iterator it = states.begin();
		while(it != states.end()){
			if(*it == s) return true;
			++it;
		}
		return false;
	};
	
	inline timed_state* get_state(int number){
		if(number < states.size())
			return states[number];
		return 0;
	};
	
	inline int get_number(timed_state* state){
		int i = 0;
		for(state_list::const_iterator it = states.begin(); it != states.end();	++it){
			if(*it == state) return i;
			i++;
		}
		return -1;
	};
	
	inline int num_states(){
		return states.size();
	};
	
	inline char get_alph_char(int i){
		return input->get_symbol(i);
	};

	inline char get_alph_int(char c){
		return input->get_int(c);
	};
	
	inline const timed_input* get_input(){
		return input;
	};
};

class timed_state{
private:
	interval_set* targets;
	
	inline void pre_split(timed_state* old_target, timed_state* new_target);
	inline void un_pre_split(timed_state* old_target);
	inline void recurse_merge(timed_state* old_target, timed_state* new_target);
	inline void recurse_un_merge(timed_state* old_target, timed_state* new_target);
	inline void recurse_split(interval* new_in, timed_state* old_target);
	inline void recurse_un_split(interval* new_in, timed_state* old_target);
	inline void recurse_test_merge(timed_state* old_target, timed_state* new_target);
	inline void recurse_test_split(timed_state* state);
	
	friend class timed_automaton;
	friend class state_statistics;
	friend double calculate_chi2_score(timed_state* old_target, timed_state* new_target);
	friend double calculate_chi2_score(timed_state* target);
	friend double calculate_chi2_score_time(timed_state* old_target, timed_state* new_target);
	friend double calculate_chi2_score_time(timed_state* target);
	friend pair<int, double> get_likelihood_ratio(timed_state* old_target, timed_state* new_target);
	friend pair<int, double> get_likelihood_ratio(timed_state* target);
	friend pair<int, double> get_likelihood_ratio_time(timed_state* old_target, timed_state* new_target);
	friend pair<int, double> get_likelihood_ratio_time(timed_state* target);
  
public:
	timed_state();
	timed_state(timed_state* state);    
	~timed_state();
	
	state_statistics* stat;

	void create_states();	

	string to_str(timed_automaton*);
	string to_str_full(timed_automaton*);

    inline timed_state* get_target(int symbol, int time) const{
    	return get_interval_from_set(targets[symbol], time)->get_target();
    };
    
    inline interval_set& get_intervals(const int symbol) const{
    	return targets[symbol];
    };
    
    inline interval* get_interval(int symbol, int time) const{
    	return get_interval_from_set(targets[symbol], time);
    };
    
    void add_tail(timed_tail* tail);
    void del_tail(timed_tail* tail);

    void point(int symbol, int time, timed_state* target);
    void undo_point(int symbol, int time, timed_state* target);
    double test_point(int symbol, int time, timed_state* target);
    
    void split(int symbol, int time);
    void undo_split(int symbol, int time);
    double test_split(int symbol, int time);
    
	void mark(interval*, timed_tail* tail);
	void un_mark(interval*, timed_tail* tail);	
	void clear_marked(interval*);
};

#endif /* TIMED_AUTOMATON_H_*/
