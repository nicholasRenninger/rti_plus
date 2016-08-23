/*
 *  RTI (real-time inference)
 *  Interval.h, the header file for the interval time structure
 *  The timing is from begin to end, inclusive
 *  The structure contains a list of tails (suffixes of timed strings)
 *  And a target state
 *
 *  Essentially it is a transition of a real-time automaton as is:
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

#ifndef _INTERVAL_H_
#define _INTERVAL_H_


#include <map>
#include "tail.h"
using namespace std;

class interval;
class timed_state;

typedef map<int, interval*> interval_set;
typedef map<int, interval*>::iterator interval_it;
typedef map<int, interval*>::reverse_iterator interval_rit;
typedef map<int, interval*>::const_iterator const_interval_it;

void split_set(interval_set&, int time);
void undo_split_set(interval_set&, int time);

void create_interval_set(interval_set& i_set);
void delete_interval_set(interval_set& i_set);
static inline interval* get_interval_from_set(const interval_set& intervals, int time)
{
	const_interval_it it = intervals.lower_bound(time);
	if(it == intervals.end()){ return (*intervals.rbegin()).second; }
	//assert((*it).second->get_begin() <= time && (*it).second->get_end() >= time);
	return (*it).second;
};

class interval{
private:
	int begin;        // inclusive
	int end;          // inclusive
	tail_set tails;   // tails
	timed_state* to;  // target state
	
	int num_marked;

	friend class timed_state;

	friend void split_set(interval_set&, int time);
	friend void undo_split_set(interval_set&, int time);

	friend interval_set& create_interval_set();
	friend interval* get_interval(const interval_set&, int time);

public:
	double probability;
	tail_set undo_tails;
	timed_state* undo_to;

	/* constructor */
	interval(int b, int e);
	
	/* get methods */
	inline timed_state* get_target() const{
		return to;
	};
	
	inline void set_target(timed_state* state){
		to = state;
	};
	
	inline int get_begin() const{
		return begin;
	};
	
	inline int get_end() const{
		return end;
	};
	
	inline tail_set& get_tails(){
		return tails;
	};
	
	inline void add_tail(timed_tail* tail){
		add_tail_to_set(tails, tail);
	};
	
	inline void del_tail(timed_tail* tail){
		del_tail_from_set(tails, tail);
	};
	
	inline bool contains_tail(timed_tail* tail){
		return contains_tail_in_set(tails, tail);
	};

	inline bool is_empty() const{
		return tails.empty();
	};
	
	inline void add_marked(){
		++num_marked;
	};
	
	inline void del_marked(){
		--num_marked;
	};
	
	inline int get_num_marked(){
		return num_marked;
	};
};

#endif /* _INTERVAL_H_ */
