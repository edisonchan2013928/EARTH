#The following lines are used for compiling our code.
g++ -c Network.cpp -w -o Network.o
g++ -c SP.cpp -w -o SP.o
g++ -c range_tree.cpp -w -o range_tree.o
g++ -c count.cpp -w -o count.o
g++ -c K_function.cpp -w -o K_function.o

g++ main.cpp -O3 -o main Network.o SP.o range_tree.o count.o K_function.o
exit #Remove it if you want to run our code.

#These are the parameters for calling our code.
#our_model.network_fileName = argv[1]; //the input data file 
#our_model.method = atoi(argv[2]); //method = 1: RQS, method = 2: SPS, method = 3: EAR 
#our_model.s_threshold = atof(argv[3]); //the spatial threshold (e.g., 1000) in terms of meters
#our_model.t_threshold = atof(argv[4]); //the temporal threshold (e.g., 3) in terms of days

dir="./Datasets/" #directory name of the input file
s_threshold=100 #spatial threshold = 100 meters
t_threshold=7 #temporal threshold = 7 days
method=2 #the EAR method
dataset="Seattle_network" #the input file name

./main $dir$dataset $method $s_threshold $t_threshold
