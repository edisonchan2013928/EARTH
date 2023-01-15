g++ -c Network.cpp -w -o Network.o
g++ -c SP.cpp -w -o SP.o
g++ -c range_tree.cpp -w -o range_tree.o
g++ -c count.cpp -w -o count.o
g++ -c K_function.cpp -w -o K_function.o

g++ main.cpp -O3 -o main Network.o SP.o range_tree.o count.o K_function.o
exit #Remove it if you want to run our code.

#These are the parameters for calling our code.
#our_model.network_fileName = argv[1]; //the input data file
#our_model.out_fileName = argv[2]; //the output data file
#our_model.method = atoi(argv[3]); //method = 1: RQS, method = 2: SPS, method = 3: EAR, method = 4: MTS, method = 5: EARTH
#our_model.L = atoi(argv[4]); //number of datasets (with L - 1 randomly generated datasets)
#our_model.M = atoi(argv[5]); //number of spatial thresholds
#our_model.T = atoi(argv[6]); //number of temporal thresholds
#our_model.start_s_tau = atof(argv[7]); #the initial spatial threshold (e.g., 200) in terms of meters
#our_model.incr_s_tau = atof(argv[8]); #the incremental spatial threshold (e.g., 200) in terms of meters
#our_model.start_t_tau = atof(argv[9]); #the initial temporal threshold (e.g., 3) in terms of days
#our_model.incr_t_tau = atof(argv[10]); #the incremental temporal threshold (e.g., 2) in terms of days

dir="./Datasets/" #directory name of the input file
out_dir="./Results/" #directory name of the output file
method=3
L=2
M=4
T=4
start_s_tau=200
incr_s_tau=200
start_t_tau=3
incr_t_tau=2
dataset="Seattle_network" #the input file name

./main $dir$dataset $out_dir$dataset"_md"$method"_L"$L"_M"$M"_T"$T $method $L $M $T $start_s_tau $incr_s_tau $start_t_tau $incr_t_tau
