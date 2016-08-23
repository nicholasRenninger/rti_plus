/*
 *  RTI (real-time inference)
 *  Timed_automaton.cpp, the source file for real-time automata, essentially just a set of states.
 *
 *  This file also contains the merge and split testing routines, which test is performed is determined by the TEST_TYPE value.
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
 
#include <gsl/gsl_cdf.h>
#include <math.h>
#include <string>
#include <stdio.h>
#include "timed_automaton.h"
#include <assert.h>

int TEST_TYPE = 0;

timed_automaton::timed_automaton(){
	root = new timed_state();
	states.push_back(root);
	input = 0;
};

timed_automaton::timed_automaton(timed_input* in){
	root = new timed_state();
	states.push_back(root);
	input = in;
	
	for(int i = 0; i < in->get_num_words(); ++i){
	 	timed_tail* tail = new timed_tail(in->get_word(i), 0, 0);
	 	timed_tail* prev_tail = tail;
	    for(int index = 1; index < in->get_word(i)->get_length(); ++index)
	    	prev_tail = new timed_tail(in->get_word(i), index, prev_tail);

		if(tail->get_symbol() != 10000) root->add_tail(tail);
	}
	root->create_states();
};

timed_automaton::~timed_automaton(){
	for(state_list::iterator it = states.begin(); it != states.end(); ++it)
		delete *it;
};

void timed_automaton::check_next_tail(interval* in, timed_tail* tail){
	assert(in->get_begin() <= tail->get_time_value());
	assert(in->get_end()   >= tail->get_time_value());
	assert(in->contains_tail(tail));
	assert(!tail->is_marked());
	if(tail->next_tail() != 0){
		assert(in->get_target() != 0);
		interval* next_in = in->get_target()->get_interval(tail->next_tail()->get_symbol(), tail->next_tail()->get_time_value());
		assert(next_in->contains_tail(tail->next_tail()));
		assert(next_in->get_begin() <= tail->next_tail()->get_time_value());
		assert(next_in->get_end()   >= tail->next_tail()->get_time_value());
		if(!contains_state(in->get_target())){
			assert(in->get_target()->stat->get_total_marks() == 0);
			assert(next_in->get_begin() == MIN_TIME);
			assert(next_in->get_end() == MAX_TIME);
			timed_automaton::check_next_tail(next_in, tail->next_tail());
		}
	} 
};

void timed_automaton::check_consistency(){
#ifdef NDEBUG
  return;
#endif
	for(state_list::iterator it1 = states.begin(); it1 != states.end(); ++it1){
		timed_state* state = *it1;
		assert(state->stat->get_total_marks() == 0);
		for(int i = 0; i < MAX_SYMBOL; ++i){
			for(interval_it it2 = state->get_intervals(i).begin(); it2 != state->get_intervals(i).end(); ++it2){
				interval* in = (*it2).second;
				assert(in->get_num_marked() == 0);
				for(const_tail_it it3 = in->get_tails().begin(); it3 != in->get_tails().end(); ++it3){
					assert(in->get_begin() <= (*it3).first);
					assert(in->get_end()   >= (*it3).first);
					assert(in->contains_tail((*it3).second));
					timed_tail* tail = (*it3).second;
					assert(!tail->is_marked());
					if(tail->next_tail() != 0){
						assert(in->get_target() != 0);
						timed_automaton::check_next_tail(in->get_target()->get_interval(tail->next_tail()->get_symbol(), tail->next_tail()->get_time_value()), tail->next_tail());
					} 
				}
			}
		}
	}
};

void timed_automaton::recursive_tree_automaton(timed_state* st, timed_state* garbage_state){
	add_state(st);
	for(int s = 0; s < MAX_SYMBOL; ++s){
		for(const_interval_it it2 = st->get_intervals(s).begin(); it2 != st->get_intervals(s).end(); ++it2){
			interval* in = (*it2).second;
			if(in->get_target() == 0) in->set_target(garbage_state);
			
			if(in->is_empty()) continue;
			recursive_tree_automaton(in->get_target(), garbage_state);
		}
	}
};

void timed_automaton::tree_automaton(){
	timed_state* garbage_state = new timed_state();
	for(int i = 0; i < MAX_SYMBOL; ++i){
		garbage_state->point(i, 0, garbage_state);
	}

	for(state_it it = get_states().begin(); it != get_states().end(); ++it){
		timed_state* st = *it;
		for(int s = 0; s < MAX_SYMBOL; ++s){
			for(const_interval_it it2 = st->get_intervals(s).begin(); it2 != st->get_intervals(s).end(); ++it2){
				interval* in = (*it2).second;
				if(in->get_target() == 0) in->set_target(garbage_state);
				else{
					if(contains_state(in->get_target()) || in->is_empty()) continue;
					recursive_tree_automaton(in->get_target(), garbage_state);
				}
			}
		}
	}
	add_state(garbage_state);
};

int timed_automaton::recursive_total_num_states(timed_state* st){
	int result = 1;
	for(int s = 0; s < MAX_SYMBOL; ++s){
		for(const_interval_it it2 = st->get_intervals(s).begin(); it2 != st->get_intervals(s).end(); ++it2){
			interval* in = (*it2).second;
			if(in->get_target() == 0 || in->get_target() == st) continue;
			result += recursive_total_num_states(in->get_target());
		}
	}
	return result;
};

int timed_automaton::total_num_states(){
	int result = 0;
	for(state_it it = get_states().begin(); it != get_states().end(); ++it){
		result++;
		timed_state* st = *it;
		for(int s = 0; s < MAX_SYMBOL; ++s){
			for(const_interval_it it2 = st->get_intervals(s).begin(); it2 != st->get_intervals(s).end(); ++it2){
				interval* in = (*it2).second;
				if(in->get_target() == 0) continue;
				if(contains_state(in->get_target())) continue;

				result += recursive_total_num_states(in->get_target());
			}
		}
	}
	return result;
};

int timed_automaton::get_size(){
	int result = 0;
	for(state_it it = get_states().begin(); it != get_states().end(); ++it){
		for(int s = 0; s < MAX_SYMBOL; ++s){
			timed_state* prev_state = 0;
			bool first = true;
			for(const_interval_it it2 = (*it)->get_intervals(s).begin(); it2 != (*it)->get_intervals(s).end(); ++it2){
				interval* in = (*it2).second;
				if(!first && in->get_target() == prev_state) continue;
				result++;
				prev_state = in->get_target();
				first = false;
			}
		}
	}
	return result;
};

void timed_automaton::garbage_automaton(){
	timed_state* root = (*states.begin());
	for(int i = 0; i < MAX_SYMBOL; ++i){
		root->point(i, 0, root);
	}
};

void timed_automaton::from_file(FILE * str){
	assert(states.size() == 1);
	for(int i = 0; i < MAX_SYMBOL; ++i){
		assert(get_number(root->get_target(i, 1)) == -1);
	}
	
	timed_state* garbage_state = new timed_state();
	for(int i = 0; i < MAX_SYMBOL; ++i){
		garbage_state->point(i, 0, garbage_state);
	}
	
	while(!feof(str)){
		int source_state, begin_time, end_time, target_state, num_strings;
		float probability;
		char symbol;
		int num_conversions = fscanf(str, "%d %c [%d, %d]->%d #%d p=%f\n", &source_state, &symbol, &begin_time, &end_time, &target_state, &num_strings, &probability);
		if(num_conversions != 7) return;
		
		if(get_alph_int(symbol) != -1){
			while(get_state(source_state) == 0) add_state(new timed_state());
			timed_state* s = get_state(source_state);
			while(target_state != -1 && get_state(target_state) == 0) add_state(new timed_state());
			timed_state* t = get_state(target_state);
			if(t == 0) t = garbage_state;
			
			if(begin_time > MAX_TIME) continue;
			if(end_time > MAX_TIME) end_time = MAX_TIME;
			
			interval* in = s->get_interval(get_alph_int(symbol), begin_time);
			if(in->get_begin() != begin_time) s->split(get_alph_int(symbol), begin_time - 1); 
			in = s->get_interval(get_alph_int(symbol), begin_time);
			if(in->get_end() != end_time) s->split(get_alph_int(symbol), end_time);
			in = s->get_interval(get_alph_int(symbol), begin_time);
			assert(in->get_begin() == begin_time);
			assert(in->get_end() == end_time);
			
			s->point(get_alph_int(symbol), begin_time, t);
			assert(s->get_target(get_alph_int(symbol), begin_time) == t);
		}
	}
	
	for(state_it it = get_states().begin(); it != get_states().end(); ++it){
		timed_state* st = *it;
		for(int s = 0; s < MAX_SYMBOL; ++s){
			for(const_interval_it it2 = st->get_intervals(s).begin(); it2 != st->get_intervals(s).end(); ++it2){
				if(get_number((*it2).second->get_target()) == -1 && (*it2).second->get_target() != garbage_state){
					st->point(s, (*it2).second->get_begin(), garbage_state);
				}
			}
		}
	}
	add_state(garbage_state);
};

string timed_automaton::to_str_full(){
	ostringstream ostr;
	for(state_list::const_iterator it = states.begin(); it != states.end(); ++it)
		ostr << (*it)->to_str_full(this);
	return ostr.str();
};

string timed_automaton::to_str(){
	ostringstream ostr;
	for(state_list::const_iterator it = states.begin(); it != states.end(); ++it)
		ostr << (*it)->to_str(this);
	return ostr.str();
};

string timed_state::to_str_full(timed_automaton* ta){
	ostringstream ostr;
	ostr << ta->get_number(this) << " prob: symbol= ";
	for(int i = 0; i < MAX_SYMBOL; ++i){
		ostr << stat->get_symbol_counts(i) << " ";
	}
	ostr << " time= ";
	for(int i = 0; i < NUM_HISTOGRAM_BARS; ++i){
		ostr << stat->get_time_counts(i) << " ";
	}
	ostr << endl;
	for(int i = 0; i < MAX_SYMBOL; ++i){
		for(const_interval_it it = targets[i].begin(); it != targets[i].end(); ++it){
			ostr << ta->get_number(this) << " "  << i
			<< " [" << (*it).second->get_begin()
			<< ", ";
			ostr << (*it).second->get_end()
				<< "]->" << ta->get_number((*it).second->get_target()) << endl;
		}
	}
	return ostr.str();
};

string timed_state::to_str(timed_automaton* ta){
	ostringstream ostr;
	double total_size = 0.0;
	for(int i = 0; i < MAX_SYMBOL; ++i){
		for(const_interval_it it = targets[i].begin(); it != targets[i].end(); ++it){
		  total_size += (*it).second->tails.size();
		}
	}
	for(int i = 0; i < MAX_SYMBOL; ++i){
		const_interval_it prev_it = targets[i].begin();
		int prev_time = -1;
		int prev_size = 0;
		for(const_interval_it it = targets[i].begin(); it != targets[i].end(); ++it){
		  if((*it).second->tails.size() != 0){
				if((*prev_it).second->get_target() != (*it).second->get_target()){
					if(prev_size != 0){
						ostr << ta->get_number(this) << " "  << ta->get_alph_char(i)
						<< " [" << (*prev_it).second->get_begin()
						<< ", ";
						ostr << prev_time
							<< "]->" << ta->get_number((*prev_it).second->get_target());
						ostr << " #" << prev_size << " p=" << ((double)prev_size / (double)total_size) << "\n";
					}
					prev_size = (*it).second->tails.size();
					prev_time = (*it).second->get_end();
					prev_it = it;
				} else {
					prev_size += (*it).second->tails.size(); 
					prev_time = (*it).second->get_end();
				}
			}
		}
		if(prev_size != 0){
			ostr << ta->get_number(this) << " "  << ta->get_alph_char(i)
			<< " [" << (*prev_it).second->get_begin()
			<< ", ";
			ostr << prev_time
				<< "]->" << ta->get_number((*prev_it).second->get_target());
			ostr << " #" << prev_size << " p=" << ((double)prev_size / (double)total_size) << "\n";
		}
	}
	return ostr.str();
};

timed_state::timed_state(){
	targets = new interval_set[MAX_SYMBOL];
	for(int i = 0; i < MAX_SYMBOL; ++i)
		create_interval_set(targets[i]);
		
	stat = new state_statistics();
}

timed_state::timed_state(timed_state* state){
	stat = new state_statistics();

	targets = new interval_set[MAX_SYMBOL];
	for(int i = 0; i < MAX_SYMBOL; ++i)
		create_interval_set(targets[i]);
	
	for(int i = 0; i < MAX_SYMBOL; ++i){
		for(interval_it it = state->get_intervals(i).begin(); it != state->get_intervals(i).end(); ++it){
			interval* in = (*it).second;
			if(in->get_end() != MAX_TIME)
				split_set(targets[i], in->get_end());
		}
	}

	for(int i = 0; i < MAX_SYMBOL; ++i){
		interval_it it1 = get_intervals(i).begin();
		interval_it it2 = state->get_intervals(i).begin();
		while(it1 != get_intervals(i).end()){
			(*it1).second->to = (*it2).second->to;
			++it1;
			++it2;
		}
	}
}

void timed_state::create_states(){
	for(int i = 0; i < MAX_SYMBOL; ++i){
		for(interval_it it = targets[i].begin(); it != targets[i].end(); ++it){
			interval* in = (*it).second;
			assert(in->get_target() == 0);
			if(!in->is_empty()){
				in->to = new timed_state();
				for(const_tail_it it = in->get_tails().begin(); it != in->get_tails().end(); ++it){
					timed_tail* tail = (*it).second;
					if(tail->next_tail() != 0) in->to->add_tail(tail->next_tail());
				}
				in->to->create_states();
			}
		}
	}	
}

timed_state::~timed_state(){
	for(int i = 0; i < MAX_SYMBOL; ++i)
		delete_interval_set(targets[i]);
	delete[] targets;

	delete stat;
};

void timed_state::add_tail(timed_tail* tail){
	get_interval_from_set(targets[tail->get_symbol()], tail->get_time_value())->add_tail(tail);
	stat->add_count(tail);
};

void timed_state::del_tail(timed_tail* tail){
   	get_interval_from_set(targets[tail->get_symbol()], tail->get_time_value())->del_tail(tail);
	stat->del_count(tail);
};

void timed_state::pre_split(timed_state* old_target, timed_state* new_target){
	for(int i = 0; i < MAX_SYMBOL; ++i){
		if((*old_target->get_intervals(i).begin()).first != MAX_TIME){
			cerr << (*old_target->get_intervals(i).begin()).first << endl;
			cerr << (*new_target->get_intervals(i).begin()).first << endl;
			cerr << (*old_target->get_intervals(i).begin()).second->is_empty() << endl;
			cerr << (*new_target->get_intervals(i).begin()).second->is_empty() << endl;
			cerr << (*old_target->get_intervals(i).begin()).second->get_target() << endl;
			cerr << (*new_target->get_intervals(i).begin()).second->get_target() << endl;
			assert(0);
		}
		for(interval_it it = new_target->get_intervals(i).begin(); it != new_target->get_intervals(i).end(); ++it){
			interval* new_in = (*it).second;
			if(new_in->get_end() != MAX_TIME)
				old_target->split(i, new_in->get_end());
		}
	}
	for(int i = 0; i < MAX_SYMBOL; ++i){
		interval_it it1 = old_target->get_intervals(i).begin();
		interval_it it2 = new_target->get_intervals(i).begin();
		while(it1 != old_target->get_intervals(i).end()){
			interval* old_in = (*it1).second;
			interval* new_in = (*it2).second;
			assert(old_in->get_end() == new_in->get_end());
			assert(old_in->get_begin() == new_in->get_begin());
			if(old_in->to != 0 && new_in->to != 0) pre_split(old_in->to, new_in->to);
			++it1;
			++it2;
		}
	}
};

void timed_state::un_pre_split(timed_state* old_target){
	for(int i = 0; i < MAX_SYMBOL; ++i){
		interval_it it1 = old_target->get_intervals(i).begin();
		while(it1 != old_target->get_intervals(i).end()){
			interval* old_in = (*it1).second;
			if(old_in->to != 0) un_pre_split(old_in->to);
			++it1;
		}
	}
	for(int i = 0; i < MAX_SYMBOL; ++i){
		interval* old_in = (*old_target->get_intervals(i).begin()).second;
		while(old_in->get_end() != MAX_TIME){
			old_target->undo_split(i, old_in->get_end());
			old_in = (*old_target->get_intervals(i).begin()).second;
		}
		assert((*old_target->get_intervals(i).begin()).first == MAX_TIME);
	}			
};

void timed_state::recurse_split(interval* new_in, timed_state* old_target){
	timed_state* new_target = new_in->get_target();
	for(const_tail_it it = new_in->get_tails().begin(); it != new_in->get_tails().end(); ++it){
		timed_tail* tail = (*it).second;
		if(tail->next_tail() != 0){
			old_target->del_tail(tail->next_tail());
			new_target->add_tail(tail->next_tail());
		}
	}
	
	for(int i = 0; i < MAX_SYMBOL; ++i){
		interval* next_new_in = new_target->get_interval(i, MAX_TIME);
		if(!next_new_in->is_empty()){
			interval* next_old_in = old_target->get_interval(i, MAX_TIME);
			if(!next_old_in->is_empty()){
				next_new_in->to = new timed_state();
				recurse_split(next_new_in, next_old_in->to);
			} else {
				next_new_in->to = next_old_in->to;
				next_old_in->to = 0;
			}			
		}
	}
};

void timed_state::recurse_un_split(interval* new_in, timed_state* old_target){
	timed_state* new_target = new_in->get_target();
	for(int i = MAX_SYMBOL - 1; i >= 0; --i){
		interval* next_new_in = new_target->get_interval(i, MAX_TIME);
		if(!next_new_in->is_empty()){
			interval* next_old_in = old_target->get_interval(i, MAX_TIME);
			if(!next_old_in->is_empty()){
				recurse_un_split(next_new_in, next_old_in->to);
				delete next_new_in->to;
			} else {
				next_old_in->to = next_new_in->to;
				next_new_in->to = 0;
			}
		}
	}
	
	for(int i = MAX_SYMBOL - 1; i >= 0; --i){
		interval* next_new_in = new_target->get_interval(i, MAX_TIME);
		interval* next_old_in = old_target->get_interval(i, MAX_TIME);
		
		next_old_in->tails.insert(next_new_in->tails.begin(), next_new_in->tails.end());
		for(tail_it it3 = next_new_in->tails.begin(); it3 != next_new_in->tails.end(); ++it3){
			old_target->stat->add_count((*it3).second);
			new_target->stat->del_count((*it3).second);
		}
	}
};

void timed_state::recurse_merge(timed_state* old_target, timed_state* new_target){
	for(int i = 0; i < MAX_SYMBOL; ++i){
		interval_it it1 = old_target->get_intervals(i).begin();
		interval_it it2 = new_target->get_intervals(i).begin();
		while(it1 != old_target->get_intervals(i).end()){
			interval* old_in = (*it1).second;
			interval* new_in = (*it2).second;
			assert(old_in->get_end() == new_in->get_end());
			assert(old_in->get_begin() == new_in->get_begin());
			
			if(!old_in->is_empty()){
				if(!new_in->is_empty()){
					recurse_merge(old_in->to, new_in->to);
				} else {
					new_in->to = old_in->to;
					old_in->to = 0;
				}
				new_in->tails.insert(old_in->tails.begin(), old_in->tails.end());
				for(tail_it it3 = old_in->tails.begin(); it3 != old_in->tails.end(); ++it3){
					new_target->stat->add_count((*it3).second);
				}
			}
			++it1;
			++it2;
		}
	}
};

void timed_state::recurse_un_merge(timed_state* old_target, timed_state* new_target){
	for(int i = MAX_SYMBOL - 1; i >= 0; --i){
		interval_rit it1 = old_target->get_intervals(i).rbegin();
		interval_rit it2 = new_target->get_intervals(i).rbegin();
		while(it1 != old_target->get_intervals(i).rend()){
			interval* old_in = (*it1).second;
			interval* new_in = (*it2).second;
			assert(old_in->get_end() == new_in->get_end());
			assert(old_in->get_begin() == new_in->get_begin());
			
			if(!old_in->is_empty()){
				for(tail_it it3 = old_in->tails.begin(); it3 != old_in->tails.end(); ++it3){
					new_in->del_tail((*it3).second);
					new_target->stat->del_count((*it3).second);
				}
				if(!new_in->is_empty()){
					recurse_un_merge(old_in->to, new_in->to);
				} else {
					old_in->to = new_in->to;
					new_in->to = 0;
				}
			}
			++it1;
			++it2;
		}
	}
/*
	for(int i = 0; i < MAX_SYMBOL; ++i){
		interval_it it1 = old_target->get_intervals(i).begin();
		interval_it it2 = new_target->get_intervals(i).begin();
		while(it1 != old_target->get_intervals(i).end()){
			interval* old_in = (*it1).second;
			interval* new_in = (*it2).second;
			assert(old_in->get_end() == new_in->get_end());
			assert(old_in->get_begin() == new_in->get_begin());
			
			if(!old_in->is_empty()){
				if(!new_in->is_empty()){
					old_in->to = new timed_state(new_in->to);
					for(const_tail_it it = old_in->get_tails().begin(); it != old_in->get_tails().end(); ++it){
						timed_tail* tail = (*it).second;
						if(tail->next_tail() != 0){
							new_in->to->del_tail(tail->next_tail());
							old_in->to->add_tail(tail->next_tail());
						}
					}
					recurse_un_merge(old_in->to, new_in->to);
				} else {
					new_in->to = 0;
				}
			} else {
				old_in->to = 0;
			}
			++it1;
			++it2;
		}
	}*/
};

void timed_state::split(int symbol, int time){
	interval* in = get_interval(symbol, time);
	split_set(targets[symbol], time);
	interval* new_in = get_interval(symbol, time);
	assert(new_in != in && new_in->get_target() == 0);
	
	if(!new_in->is_empty()){
		if(!in->is_empty()){
			new_in->to = new timed_state();
			recurse_split(new_in, in->get_target());
		} else {
			new_in->to = in->to;
			in->to = 0;
		}
	}		
};

void timed_state::undo_split(int symbol, int time){
	interval* in = get_interval(symbol, time + 1);
	interval* new_in = get_interval(symbol, time);
	
	if(!new_in->is_empty()){
		if(!in->is_empty()){
			recurse_un_split(new_in, in->get_target());
			delete new_in->to;
		} else {
			in->to = new_in->to;
			new_in->to = 0;
		}
	}
	undo_split_set(targets[symbol], time);
};

void timed_state::point(int symbol, int time, timed_state* new_target){
	interval* in = get_interval(symbol, time);
	in->undo_tails = tail_set(in->tails);
	timed_state* old_target = in->get_target();
	assert(old_target != new_target);
	in->to = new_target;
	if(old_target != 0){
		pre_split(old_target, new_target);
		recurse_merge(old_target, new_target);
		in->undo_to = old_target;
	}
};

void timed_state::undo_point(int symbol, int time, timed_state* new_target){
	interval* in = get_interval(symbol, time);
	assert(in->to == new_target);
	timed_state* old_target = in->undo_to;
	in->undo_to = 0;
/*	for(const_tail_it it = in->undo_tails.begin(); it != in->undo_tails.end(); ++it){
		timed_tail* tail = (*it).second;
		if(tail->next_tail() != 0){
			new_target->del_tail(tail->next_tail());
			old_target->add_tail(tail->next_tail());
		}
	}*/
	if(old_target != 0){
		recurse_un_merge(old_target, new_target);
		un_pre_split(old_target);
		in->to = old_target;
	}
};

void timed_state::recurse_test_merge(timed_state* old_target, timed_state* new_target){
	if(old_target == 0 || new_target == 0) return;

	if(TEST_TYPE == 2) {
		calculate_chi2_score(old_target, new_target);
		calculate_chi2_score_time(old_target, new_target);
	} else {
		get_likelihood_ratio(old_target, new_target);
		get_likelihood_ratio_time(old_target, new_target);
	}
	
	for(int i = 0; i < MAX_SYMBOL; ++i){
		interval_it it_1 = old_target->get_intervals(i).begin();
		interval_it it_2 = new_target->get_intervals(i).begin();
		while(it_1 != old_target->get_intervals(i).end()){
			interval* old_in = (*it_1).second;
			interval* new_in = (*it_2).second;
			++it_1;
			++it_2;

			if(old_in->get_tails().size() < MIN_DATA || new_in->get_tails().size() < MIN_DATA) continue;

			recurse_test_merge(old_in->to, new_in->to);
		}
	}
};

void timed_state::recurse_test_split(timed_state* state){
	if(state == 0) return;
	
	if(TEST_TYPE == 2) {
		calculate_chi2_score(state);
		calculate_chi2_score_time(state);
	} else {
		get_likelihood_ratio(state);
		get_likelihood_ratio_time(state);
	}
	
	for(int i = 0; i < MAX_SYMBOL; ++i){
		for(interval_it it = state->get_intervals(i).begin(); it != state->get_intervals(i).end(); ++it){
			interval* in = (*it).second;
			
			if(in->get_tails().size() - in->get_num_marked() < MIN_DATA || in->get_num_marked() < MIN_DATA) continue;
			
			recurse_test_split(in->to);
		}
	}
};

double timed_state::test_point(int symbol, int time, timed_state* new_target){
	timed_state* old_target = get_interval(symbol, time)->get_target();
	if(old_target == 0) return 0.0;
	
	assert(old_target != new_target);
	
	if(TEST_TYPE == 2) initialize_consensus_test();
	else initialize_likelihood_test();
	
	get_interval(symbol, time)->to = new_target;
	pre_split(old_target, new_target);
	recurse_test_merge(old_target, new_target);
	un_pre_split(old_target);
	get_interval(symbol, time)->to = old_target;

	double p_value = 0.0;
	if(TEST_TYPE == 2) p_value = calculate_consensus_test();
	else p_value = calculate_likelihood_test();
	
	return p_value;
};

void timed_state::mark(interval* in, timed_tail* tail){
	if(tail->is_marked()) return;
	stat->mark(tail);
	in->add_marked();
	tail->mark();
	if(tail->next_tail() != 0)
		in->to->mark(in->to->get_interval(tail->next_tail()->get_symbol(), tail->next_tail()->get_time_value()), tail->next_tail());
};

void timed_state::un_mark(interval* in, timed_tail* tail){
	if(!tail->is_marked()) return;
	stat->unmark(tail);
	in->del_marked();
	tail->un_mark();
	if(tail->next_tail() != 0)
		in->to->un_mark(in->to->get_interval(tail->next_tail()->get_symbol(), tail->next_tail()->get_time_value()), tail->next_tail());
};

void timed_state::clear_marked(interval* in){
	for(const_tail_it it = in->tails.begin(); it != in->tails.end(); ++it)
		un_mark(in, (*it).second);
};

double timed_state::test_split(int symbol, int time){
	if(TEST_TYPE == 2) initialize_consensus_test();
	else initialize_likelihood_test();

	interval* in = get_interval(symbol, time);
	timed_state* target = in->get_target();
	if(target == 0) return 0.0;
	for(const_tail_it it = in->get_tails().begin(); it != in->get_tails().end(); ++it){
		if((*it).second->get_time_value() <= time)
			mark(in, (*it).second);
		else
			assert(!(*it).second->is_marked());
	}

	recurse_test_split(target);
	
	double p_value = 0.0;
	if(TEST_TYPE == 2) p_value = calculate_consensus_test();
	else p_value = calculate_likelihood_test();

	return p_value;
};

