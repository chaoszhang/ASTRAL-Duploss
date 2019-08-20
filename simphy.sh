rm duploss/*/*
#simphy -sl f:25 -rs 20 -rl f:20 -rg 1 -sb f:0.000000005 -sd f:0  -st ln:21.1,0.5  -so f:1 -si f:1 -sp  f:2000000 -su ln:-17.27461,0.6931472 -hh f:1 -hs ln:1.5,1 -hl ln:1.551533,0.6931472 -hg ln:1.5,1 -cs 9644 -v 3 -o ASTRALIII -ot 0 -op 1 -lb f:0.0000000002 -ld f:0.0000000002

simphy -sl f:25 -rs 1 -rl f:5000 -rg 1 -sb f:0.000000005 -sd f:0  -st ln:21.1,0.5  -so f:1 -si f:1 -sp  f:200000000 -su ln:-17.27461,0.6931472 -hh f:1 -hs ln:1.5,1 -hl ln:1.551533,0.6931472 -hg ln:1.5,1 -cs 9644 -v 3 -o duploss -ot 0 -op 1 -lb f:0.000000002 -ld f:0.000000002 -lt f:0.000000002

cat duploss/1/g_trees* > genetrees.tre
python preprocessing.py genetrees.tre > genetrees.txt
g++ -std=c++11 GenetreeAnnotator.cpp
./a.out
java -D"java.library.path=/" -jar ASTRAL-MP/astral.5.14.2_original.jar -i breaked_trees.tre -o ASTRAL-MP/output_break.tre
java -D"java.library.path=ASTRAL-MP/lib" -jar ASTRAL-MP/astral.5.14.2.jar -i cluster_trees.tre -o ASTRAL-MP/output.tre