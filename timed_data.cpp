/*
 *  RTI (real-time inference)
 *  Timed_data.cpp, the source file for timed data, including reading from file.
 *
 *	The file format is:
 *  int int			(number_of_strings size_of_alphabet)
 *  int char int char int char int ... char int				(length_of_string symbol1 time_delay1 s2 t2 .. sn tn)
 *
 *  All ints are positive integers, greater or equal to 0
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

#include "timed_data.h"
#include <assert.h>
#include <stdlib.h>
#include <set>

int MAX_SYMBOL = 2;
int MIN_TIME = 0;
int MAX_TIME = 1000000;
int NUM_WORDS = 0;
int TOTAL_NUM_SYMBOLS = 0;

int NUM_HISTOGRAM_BARS = 4;
int TIME_IQR25 = 0;
int TIME_IQR50 = 0;
int TIME_IQR75 = 0;

timed_input::timed_input(istream &str){
	str >> num_words >> alph_size;
	set<int> time_points;
	NUM_WORDS = num_words;
	MAX_SYMBOL = alph_size;
	int time_sup = -1;
	alphabet = new char[alph_size];
	words = new timed_word*[num_words];
	int current_alphabet_size = 0;
	for(int line = 0; line < num_words; ++line)
	{
		timed_word* word = new timed_word();
	    str >> word->length;
	    word->symbols       = (int*)malloc((word->length + 1)*sizeof(int));
	    word->char_symbols  = (char*)malloc((word->length + 1)*sizeof(char));
	    word->time_values   = (int*)malloc((word->length + 1)*sizeof(int));
	    int index;
	    int time_sum = 0;
	    for(index = 0; index < word->length; ++index)
	    {
			TOTAL_NUM_SYMBOLS++;
			char c;
			str >> c;
			word->char_symbols[index] = c;
			str >> word->time_values[index];
			time_points.insert(word->time_values[index]);
			time_sum += word->time_values[index];
			bool found = 0;
			for(unsigned i = 0; i < current_alphabet_size; ++i)
			{
				if(alphabet[i] == c)
				{
					word->symbols[index] = i;
					found = 1;
					break;
				}
			}
			if(found == 0)
			{
				alphabet[current_alphabet_size] = c;
				word->symbols[index] = current_alphabet_size;
				current_alphabet_size++;
				assert(current_alphabet_size <= MAX_SYMBOL);
			}
	    }
	    word->symbols[index] = num_words;
	    word->char_symbols[index] = '\0';
	    word->time_values[index] = time_sum;
	    words[line] = word;
	}
	int number = 0;
	int iq25 = (int)(((double)time_points.size()) / 4.0);
	int iq50 = (int)(((double)time_points.size()) / 2.0);
	int iq75 = (int)(((double)time_points.size() * 3.0) / 4.0);
	for(set<int>::iterator it = time_points.begin(); it != time_points.end(); ++it){
		if(number == iq25) TIME_IQR25 = *it;
		if(number == iq50) TIME_IQR50 = *it;
		if(number == iq75) TIME_IQR75 = *it;
		if(*it > time_sup) time_sup = *it;
		number++;
	}
	MAX_TIME = time_sup;
};

timed_input::~timed_input(){
	for(int line = 0; line < num_words; ++line)
		delete words[line];
	delete[] words;
	delete[] alphabet;
};

const string timed_input::to_str() const{
  ostringstream ostr;
	ostr << num_words  << " " << alph_size << "\n";
	for(int line = 0; line < num_words; ++line)
	{
	    timed_word* word = words[line];
	    ostr << word->length << " ";
	    for(int index = 0; index < word->length; ++index)
	    {
	      ostr << word->char_symbols[index] << " ";
	      ostr << word->time_values[index] << " ";
			}
			ostr << "\n";
	  }
    return ostr.str();
};

timed_word::timed_word(){
	symbols = 0;
	time_values = 0;
	char_symbols = 0;
	length = 0;
};

