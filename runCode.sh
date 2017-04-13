g++ -o output -fopenmp main.cpp
./output | tee /tmp/tmp; less /tmp/tmp | grep "time wa"
echo "Running on "  $OMP_NUM_THREADS "threads."
./display_terrain.py
