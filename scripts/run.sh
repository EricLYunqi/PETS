#! /usr/bin

dataname=$1
graphtype=$2

RED='\033[0;31m'

compile_file="build/Release/compile.log"
output_file="output/${dataname}.log"

cmake -B build/Release -DCMAKE_BUILD_TYPE=Release

nohup cmake --build build/Release --clean-first -- -j $(nproc) > "${compile_file}">&1 &

oldpid=$!
wait $oldpid
oldpid=$?
echo "cmake exit with $oldpid"

if [ $oldpid -eq 0 ]
then
    nohup ./build/Release/pets ${dataname} ${graphtype} > "${output_file}" 2>&1 &
    echo $!
    wait $!
    echo $? >> ${output_file}
    echo "aha"
else
    echo -e "${RED}Fail to compile, check ${compile_file}"
fi