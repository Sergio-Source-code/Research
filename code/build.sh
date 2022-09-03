if [ ! -d "bin" ]; then
    mkdir "bin"
fi
cd bin
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 4
cd ../

#if [ ! -d "Debug" ]; then
#    mkdir "Debug"
#fi
#cd Debug
#cmake .. -DCMAKE_BUILD_TYPE=Debug
#make -j 4

