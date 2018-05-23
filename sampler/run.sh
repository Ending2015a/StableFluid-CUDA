nx=500
ny=500
dx=0.1
dy=0.1
rate=3
array=(350  350 490  430\
    10 350 150 430\
    180 400 320 480)
output=map.sr

./main $nx $ny $dx $dy $rate $output ${array[@]}
