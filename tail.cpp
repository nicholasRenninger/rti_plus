/*
 *  RTI (real-time inference)
 *  Tail.cpp, the source file for timed tails, or timed string suffixes
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
 
#include "tail.h"
#include <assert.h>

void add_tail_to_set(tail_set& tails, timed_tail* tail){
	tails.insert(pair<int, timed_tail*>(tail->get_time_value(), tail));
};

void del_tail_from_set(tail_set& tails, timed_tail* tail){
	pair<tail_it, tail_it> it_pair = tails.equal_range(tail->get_time_value());
	for(tail_it it = it_pair.first; it != it_pair.second; ++it){
		if((*it).second == tail){
			tails.erase(it);
			return;
		}
	}
	assert(0);
};

bool contains_tail_in_set(tail_set& tails, timed_tail* tail){
	pair<tail_it, tail_it> it_pair = tails.equal_range(tail->get_time_value());
	for(tail_it it = it_pair.first; it != it_pair.second; ++it){
		if((*it).second == tail){
			return true;
		}
	}
	return false;
};

timed_tail::timed_tail(timed_word* w, int i, timed_tail* t){
	word = w;
	index = i;
	length = w->get_length() - index;
	next = 0;
	prev = t;
	if(t != 0) t->next = this;
	marker = false;
};

timed_tail::~timed_tail(){
	if(next != 0)
		delete next;
};

const string timed_tail::to_str() const{
	ostringstream ostr;
	ostr << "(" << get_symbol() << "," << get_time_value() << ")";
	const timed_tail *nt = next;
	while(nt != 0){
		ostr << nt->get_symbol();
		nt = nt->next;
	}
	return ostr.str();
};
