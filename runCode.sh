g++ -o output -fopenmp main.cpp
./output | tee /tmp/tmp; less /tmp/tmp | grep "time wa"
echo $OMP_NUM_THREADS
sed 's/,$//g' path.txt > /tmp/path; mv /tmp/path path.txt
chmod 777 display_terrain.py
./display_terrain.py
~                      
