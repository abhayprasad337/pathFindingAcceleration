g++ -o output -fopenmp main.cpp
./output | tee /tmp/tmp; less /tmp/tmp | grep "time wa"
echo $OMP_NUM_THREADS
./display_terrain.py
