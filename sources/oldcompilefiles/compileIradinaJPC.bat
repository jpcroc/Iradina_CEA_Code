@echo off
echo %0 in current dir=%cd%
echo ...
set HOME=%userprofile%
set IRADINA_SRC=%cd%\iradina_cea

cd %IRADINA_SRC%

echo try: set IRA
echo try: where make
echo try: where git
echo try: where gcc

set IRA

set PATH=^
C:\windows\system32;^
C:\windows;^
C:\Program Files\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\bin;^
C:\Program Files (x86)\CMake 2.8\bin;^
C:\Users\jc148490\AppData\Local\Programs\Git\cmd;^
C:\Users\jc148490\AppData\Local\Programs\Git\mingw64\bin;^
C:\Users\jc148490\AppData\Local\Programs\Git\bin;^
C:\Users\jc148490\AppData\Local\Programs\Git\usr\bin;^
C:\Program Files (x86)\Subversion\bin;^

echo *
echo ********************************************
echo * compile iradina
echo ********************************************
echo *
@echo on

where make
where uname

make -f makefile_cea clean
make -f makefile_cea iradina
make -f makefile_cea install
make -f makefile_cea installGUI

