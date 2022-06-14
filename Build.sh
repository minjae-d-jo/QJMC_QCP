#set -e

srcDir=Src
libDir=Lib
binDir=Bin


function build {
	echo -e ":: $1"
	g++ -std=c++11 -Wall \
		-O3 -flto\
		-I $srcDir -I $libDir \
		-o $binDir/$1 \
		$srcDir/$1.cpp
}

function testBuild {
	echo -e ":: $1"
	g++ -std=c++11 -Wall \
		-g -fsanitize=address -fno-omit-frame-pointer \
		-lboost_unit_test_framework -lboost_iostreams \
		-I $srcDir -I $libDir \
		-o $binDir/$1 \
		$srcDir/$1.cpp
}

function debugBuild {
	echo -e ":: $1"
	g++ -DKCORE_DEBUG \
		-std=c++11 -Wall \
		-g -fsanitize=address -fno-omit-frame-pointer \
		-lboost_unit_test_framework -lboost_iostreams \
		-I $srcDir -I $libDir \
		-o $binDir/$1 \
		$srcDir/$1.cpp
}


build QJMC_ManyBodySteady_8
