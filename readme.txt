
RTI - real-time inference, code for the algorithm from: 

Sicco Verwer and Mathijs de Weerdt and Cees Witteveen (2007),
An algorithm for learning real-time automata,
In Maarten van Someren and Sophia Katrenko and Pieter Adriaans (Eds.),
Proceedings of the Sixteenth Annual Machine Learning Conference of Belgium and the Netherlands (Benelearn),
pp. 128-135.

You can use it to learn real-time automata unsupervised from unlabeled data. It is based on statistical methods and state-merging transition-splitting.

It currently runs a best first beam-like search, similar to the ed-beam algorithm for untimed automata, it can also be run using only a single greedy path by replacing he bestfirst() call by greedy() in the main routine (in searcher.cpp). It can aso be run in tet mode, in this case you can choose which merge or split to perform by yourself (replace bestfirst() with split()).

Make using:

make
make debug, includes debug information for gdb

The code requires GSL to be installed:

http://www.gnu.org/software/gsl/

Change the directories of gsl in the Makefile.

Run using:

./rti 1 0.05 filename

1 specifies the used method (1 for likelihood ratio, 2 for chi-squared)
0.05 is the significance level used in the tests
filename is a file in the following format:

int int                                     (number_of_strings size_of_alphabet)
int char int char int char int ... char int (length_of_string symbol1 time_delay1 s2 t2 .. sn tn)

see test.data for an example
test.aut us the real-time automaton used to generate this data
test.test_set is another (larger) data set generated from this automaton

More info: siccoverwer@gmail.com

Updates/code cleaning of this program will come soon.

I also have some evaluation routines for resulting real-time automata, if you need them, please ask, but it is not very clear how to evaluate these automata precisely (I use the AIC/BIC or Perplexity measures, but perhaps there are other better ones).

The output is a list of transitions of this form:
state_nr symbol [lower_time_bound, upper_time_bound]->state_nr #nr_strings p=probability_of_transition
for example:
0 0 [0, 950]->1 #1330 p=0.665
denotes a transition from state 0 (the start state) to state 1, with label 0, and clock guard [0, 950],
1330 timed strings fire this transition, which is 66,5% of all strings that reach the start state.


