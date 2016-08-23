/*
 *  RTI (real-time inference)
 *  Searcher.cpp, the header file for the search routines
 *  Currently, only a simple greedy (best-first) routine is implemented, search routines will be added later.
 *
 *  A refinement is either a point (merge), split, or color (adding of a new state) in the current real-time automaton, as in:
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

#ifndef _SEARCHER_H_
#define _SEARCHER_H_

using namespace std;

#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <map>

#include "timed_automaton.h"

extern int NODES;

class refinement;

typedef list<refinement> refinement_list;
typedef multimap<double, refinement, greater<double> > refinement_set;

timed_automaton* TA;

class refinement{
	int state;
	int target;
	int symbol;
	int time;
	
public:
	int ref_count;
	
	refinement(int s, int t, int sy, int ti);

	inline void print() const{
		if(target > -1)
			cerr << "point( " << state << " [" << symbol << ", " << time << "]->" << target << " )";
		else if(target == -1)
			cerr << " split( " << state << " [" << symbol << " , " << time << "] )";
		else
			cerr << "new( " << state << " [" << symbol << ", " << time << "]-> new )";
		cerr << endl;
	};
	
	void refine(){
		//cerr << "do : "; print();
		if(target > -1)
			TA->get_state(state)->point(symbol, time, TA->get_state(target));
		else if(target == -1)
			TA->get_state(state)->split(symbol, time);
		else
			TA->add_state(TA->get_state(state)->get_target(symbol, time));
	};
	
	void undo_refine(){
		//cerr << "undo : "; print();
		if(target > -1)
			TA->get_state(state)->undo_point(symbol, time, TA->get_state(target));
		else if(target == -1)
			TA->get_state(state)->undo_split(symbol, time);
		else
			TA->del_state(TA->get_state(state)->get_target(symbol, time));
	};
};

#endif /* _SEARCHER_H_ */
