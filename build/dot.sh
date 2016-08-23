# original by Sicco
# cat a.txt  | tr  -c "0-9. \n\r" " " | awk 'BEGIN {print "digraph a {"; }; {print "a"$1"->a"$5" [ label=\""$2"\"];" } END {print "}"}' > a.dot
# also prints time constraints and probability by Nino
cat $1  | tr  -c "0-9. \n\r" " " | awk 'BEGIN {print "digraph a {"; }; { printf("a"$1"->a"$5" [ label=\"["$3", "$4"] "$2" %.3f\"];\n", $7) } END {print "}"}' > $2
# prints only time constraints and symbol
# cat $1  | tr  -c "0-9. \n\r" " " | awk 'BEGIN {print "digraph a {"; }; { print "a"$1"->a"$5" [ label=\"["$3", "$4"] "$2"\"];\n" } END {print "}"}' > $2
