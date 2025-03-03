#! /usr/bin

dataname=$1
graphtype=$2

RED='\033[0;31m'

compile_file="build/Debug/compile.log"
output_file="output/${dataname}.log"

cmake -B build/Debug -DCMAKE_BUILD_TYPE=Debug

nohup cmake --build build/Debug --clean-first -- -j $(nproc) > "${compile_file}">&1 &

oldpid=$!
wait $oldpid
oldpid=$?
echo "cmake exit with $oldpid"

if [ $oldpid -eq 0 ]
then
    gdb --silent ./build/Debug/pets
else
    echo -e "${RED}Fail to compile, check ${compile_file}"
fi