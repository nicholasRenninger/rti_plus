/*
 *  RTI (real-time inference)
 *  Timed_data.h, the header file for timed data, including reading from file.
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
 
#ifndef _TIMED_DATA_H_
#define _TIMED_DATA_H_


class timed_input;
class timed_word;

extern int MAX_SYMBOL;
extern int MIN_TIME;
extern int MAX_TIME;
extern int NUM_HISTOGRAM_BARS;
extern int NUM_WORDS;
extern int TOTAL_NUM_SYMBOLS;

extern int TIME_IQR25;
extern int TIME_IQR50;
extern int TIME_IQR75;

#include <istream>
#include <sstream>
#include <iostream>
#include <string>
using namespace std;

class timed_input{
	char* alphabet;
	timed_word** words;
	int num_words;
	int alph_size;
	
public:
	timed_input();
	timed_input(istream& str);
	~timed_input();

	const string to_str() const;

	/* get methods */
	inline char get_symbol(int i) const{
		if(i < alph_size)
			return alphabet[i];
		return '\0';
	};
	
	inline int get_int(char c) const{
		for(int i = 0; i < alph_size; ++i)
			if(alphabet[i] == c) return i;
		return -1;
	};

	inline timed_word* get_word(int i) const{
		return words[i];
	};
	
	inline int get_num_words() const{
		return num_words;
	};

	inline int get_alph_size() const{
		return alph_size;
	};
};

class timed_word{
	int*	symbols;
	int*	time_values;
	char*	char_symbols;
	int		length;
	double	probability;

	friend class timed_input;

public:
	timed_word();

	inline const int* get_symbols() const{
		return symbols;
	};
	
	inline const char* get_char_symbols() const{
		return char_symbols;
	};

	inline const int* get_time_values() const{
		return time_values;
	};
	
	inline const int get_length() const{
		return length;
	};
	
	inline const double get_probability(){
		return probability;
	};
	
	inline void set_probability(double p){
		probability = p;
	};
};


#endif /* _TIMED_DATA_H_*/
