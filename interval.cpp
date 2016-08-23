/*
 *  RTI (real-time inference)
 *  Interval.cpp, the source file for the interval time structure
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

#include <assert.h> 
#include "interval.h"

void split_set(interval_set& intervals, int time){
	interval_it it = intervals.lower_bound(time);
	assert(it != intervals.end());
	assert((*it).second->get_begin() <= time && (*it).second->get_end() > time);

	interval* in = (*it).second;
	interval* new_in = new interval(in->get_begin(), time);
	
	
	tail_it it2 = in->tails.upper_bound(time);
	new_in->tails.insert(in->tails.begin(), it2);
	in->tails.erase(in->tails.begin(), it2);
	in->begin = time + 1;
	intervals.insert(it, pair<int, interval*>(time, new_in));
};

void undo_split_set(interval_set& intervals, int time){
	interval_it it = intervals.find(time);
	assert(it != intervals.end());

	interval* old_in = (*it).second;
	++it;
	assert(it != intervals.end());
	
	interval* in = (*it).second;
	--it;
			
	in->tails.insert(old_in->tails.begin(), old_in->tails.end());
	in->begin = old_in->get_begin();
	intervals.erase(it);
	
	delete old_in;
};

void create_interval_set(interval_set& i_set){
	interval *in = new interval(MIN_TIME, MAX_TIME);
	i_set[in->get_end()] = in;
};

void delete_interval_set(interval_set& i_set){
	for(interval_it it = i_set.begin(); it != i_set.end(); ++it)
		delete (*it).second;
	i_set.clear();
};


interval::interval(int b, int e){
	begin = b;
	end = e;
	to = 0;
	num_marked = 0;
	
	undo_to = 0;
};
