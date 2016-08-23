/*
 *  RTI (real-time inference)
 *  Tail.h, the header file for timed tails, or timed string suffixes
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
 
#ifndef _TAIL_H_
#define _TAIL_H_


#include <map>
#include <sstream>
#include "timed_data.h"

using namespace std;

class timed_tail;
typedef multimap<int, timed_tail*> tail_set;
typedef multimap<int, timed_tail*>::iterator tail_it;
typedef multimap<int, timed_tail*>::const_iterator const_tail_it;

void add_tail_to_set(tail_set&, timed_tail*);
void del_tail_from_set(tail_set&, timed_tail*);
bool contains_tail_in_set(tail_set& tails, timed_tail* tail);

class timed_tail{
	timed_word  *word;
	int          index;
	int 		 length;
	
	timed_tail *next;
	timed_tail *prev;
	
	friend class timed_input;
	
	friend void add_tail(tail_set&, timed_tail*);
	friend void del_tail(tail_set&, timed_tail*);
	
public:
	bool marker;

	/* constructs a tail from the specified arguments */
	timed_tail(timed_word *, int, timed_tail *);
	~timed_tail();

	/* used for printing */
	const string to_str() const;
	
	/* get methods */
	inline const timed_word	*get_word() const{
		return word;
	};
	
	inline int get_index() const{
		return index;
	};

	inline int get_length() const{
		return length;
	};

	inline int get_symbol() const{
		return word->get_symbols()[index];
	};
	
	inline char get_char_symbol() const{
		return word->get_char_symbols()[index];
	};

	inline const int* get_symbols() const{
		return &(word->get_symbols()[index]);
	};
	
	inline const char* get_char_symbols() const{
		return &(word->get_char_symbols()[index]);
	};

	inline int get_time_value() const{
		return word->get_time_values()[index];
	};
	
    inline int get_time_value(int i) const{
        if(i < length)
            return word->get_time_values()[index + i];
        return -1;
    };

	inline const int* get_time_values() const{
		return &(word->get_time_values()[index]);
	};
	
	inline timed_tail *next_tail() const{
		return next;
	};
	
	inline timed_tail *prev_tail() const{
		return prev;
	};
	
	inline void mark(){
		marker = true;
	};
	
	inline void un_mark(){
		marker = false;
	};
	
	inline bool is_marked(){
		return marker;
	};
};

#endif /* _TAIL_H_*/
