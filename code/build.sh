if [ ! -d "bin" ]; then
    mkdir "bin"
fi
cd bin
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 4
cd ../
