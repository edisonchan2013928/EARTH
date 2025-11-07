# EARTH
This code is the implementation of our paper "EARTH: Accelerating Spatiotemporal Network K-function-based Analytics", which is under submission to a venue.

<!-- Note that there are two folders in this Github link, which are (1) Earth and (2) Earth_plot. In the "Earth" folder, this is the implementation of all experiments for computing a spatiotemporal network K-function (Section V-B). In the "Earth_plot" folder, this is the implementation of all experiments for generating a spatiotemporal network K-function plot (Section V-C). To compile and use our code, you can refer to the "call_Earth.sh" and "call_Earth_plot.sh" files in the "Earth" and "Earth_plot" folders, respectively. These two shell script files show the details of compilation and the descriptions of different parameters. -->

# How to Compute a Spatiotemporal Network K-function?
In the "Earth" folder, this is the implementation of all experiments for computing a spatiotemporal network K-function (Section V-B). The following script shows how to compile this code.
```
g++ -c Network.cpp -w -o Network.o
g++ -c SP.cpp -w -o SP.o
g++ -c range_tree.cpp -w -o range_tree.o
g++ -c count.cpp -w -o count.o
g++ -c K_function.cpp -w -o K_function.o

g++ main.cpp -O3 -o main Network.o SP.o range_tree.o count.o K_function.o
```

After you have successfully compiled the code, you need to provide the correct parameters for calling our C++ code. 
```
our_model.network_fileName = argv[1];
our_model.method = atoi(argv[2]);
our_model.s_threshold = atof(argv[3]);
our_model.t_threshold = atof(argv[4]);
```

In the following, we discuss each parameter.

**our_model.network_fileName:** the input data file 

**our_model.method:** method = 1: RQS, method = 2: SPS, method = 3: EAR 

**our_model.s_threshold:** the spatial threshold (e.g., 1000) in terms of meters

**our_model.t_threshold:** the temporal threshold (e.g., 3) in terms of days

Here, we provide an example for running the code.
```
dir="./Datasets/" #directory name of the input file
s_threshold=100 #spatial threshold = 100 meters
t_threshold=7 #temporal threshold = 7 days
method=2 #the EAR method
dataset="Seattle_network" #the input file name

./main $dir$dataset $method $s_threshold $t_threshold
```

Please refer to the shell script file "call_Earth.sh" in the "Earth" folder for more details.

# How to Generate a Spatiotemporal Network K-function Plot?
In the "Earth_plot" folder, this is the implementation of all experiments for generating a spatiotemporal network K-function plot (Section V-C). The following script shows how to compile this code.
```
g++ -c Network.cpp -w -o Network.o
g++ -c SP.cpp -w -o SP.o
g++ -c range_tree.cpp -w -o range_tree.o
g++ -c count.cpp -w -o count.o
g++ -c K_function.cpp -w -o K_function.o

g++ main.cpp -O3 -o main Network.o SP.o range_tree.o count.o K_function.o
```
After you have successfully compiled the code, you need to provide the correct parameters for calling our C++ code. 
```
our_model.network_fileName = argv[1];
our_model.out_fileName = argv[2];
our_model.method = atoi(argv[3]);
our_model.L = atoi(argv[4]);
our_model.M = atoi(argv[5]);
our_model.T = atoi(argv[6]);
our_model.start_s_tau = atof(argv[7]);
our_model.incr_s_tau = atof(argv[8]);
our_model.start_t_tau = atof(argv[9]);
our_model.incr_t_tau = atof(argv[10]);
```

In the following, we discuss each parameter.

**our_model.network_fileName:** the input data file

**our_model.out_fileName:** the output data file

**our_model.method:** method = 1: RQS, method = 2: SPS, method = 3: EAR, method = 4: MTS, method = 5: EARTH

**our_model.L:** number of datasets (with L - 1 randomly generated datasets)

**our_model.M:** number of spatial thresholds

**our_model.T:** number of temporal thresholds

**our_model.start_s_tau:** the initial spatial threshold (e.g., 200) in terms of meters

**our_model.incr_s_tau:** the incremental spatial threshold (e.g., 200) in terms of meters

**our_model.start_t_tau:** the initial temporal threshold (e.g., 3) in terms of days

**our_model.incr_t_tau:** the incremental temporal threshold (e.g., 2) in terms of days

Here, we provide an example for running the code.
```
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
```

Please refer to the shell script file "call_Earth_plot.sh" in the "Earth_plot" folder for more details.

# Datasets
In Table II, we have provided the description of each dataset. The links of these datasets are shown as follows.

Seattle: https://data.seattle.gov/Public-Safety/SPD-Crime-Data-2008-Present/tazs-3rd5

London: https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data

New York: https://data.cityofnewyork.us/Public-Safety/Motor-Vehicle-Collisions-Crashes/h9gi-nx95

Los Angeles: https://data.lacity.org/A-Safe-City/Crime-Data-from-2010-to-2019/63jg-8b9z

In the "file_structure_description.txt" file, it provides the detailed description for storing the road network and those location data points as an input file (i.e., the input to the variable "our_model.network_fileName"). We also upload the "Testing_network" as an example to call our code.

<!-- Due to space limitations in this Github link, we do not put the datasets here. However, all these datasets are open to public, which can be downloaded from the webpages (see the references in the "Ref." column of Table II). In the "file_structure_description.txt" file, it provides the detailed description for storing the road network and those location data points as an input file (i.e., the input to the variable "our_model.network_fileName"). We also upload the "Testing_network" as an example to call our code.-->
