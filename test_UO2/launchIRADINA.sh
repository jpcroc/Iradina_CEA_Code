set -x
studyDir=/home/croc/BCA/GUIwrk/IRADINAGUI_WORKDIR/test_UO2
iradinaExe=/home/croc/BCA/GUIwrk/iradinaGUI/iradinaCode/iradina_linux64.exe
cd $studyDir
tree
ls -alt
( $iradinaExe -p 9 -data ../data -c ./Configuration.in  | tee ./iradina.log )&
# to see iradina progress
# tail -f /home/croc/BCA/GUIwrk/IRADINAGUI_WORKDIR/test_UO2/iradina.log
