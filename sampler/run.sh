nx=500
ny=500
dx=0.01
dy=0.01
rate=3
array=(350  350 490  490\
    10 350 150 490)
output=map.sr

./main $nx $ny $dx $dy $rate $output ${array[@]}
