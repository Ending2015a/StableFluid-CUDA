#!/usr/bin/python
import os

simulate = False
render = True
generate = True

input = 'map.sr'
folder = 'output'
output = 'result'
avi = 'result_3'
steps = '1000'
dt = '0.03'
nu = '1'
rho = '1000'
max_iter = '1000'
tol = '1e-7'
theme = 'Nightmare'

width = '7680'
height = '4320'

#width = 1920
#height = 1080

resize_width = 3840
resize_height = 2160



############

if simulate:

    os.system("rm ./{}/*".format(folder))
    print("simulating ...")

    output_path = os.path.join(folder, output)
    os.system("time ./main {} {} {} {} {} {} {} {}".format(input, output_path, steps, dt, nu, rho, max_iter, tol))


###############

if render:
    print("rendering ...")

    mapsr = [ os.path.join(folder, x) for x in os.listdir(folder) if x.endswith('.sr')]

    mapsr = sorted(mapsr)
    for idx, sr in enumerate(mapsr):
        print('Rendering {} frame'.format(idx))
        os.system('./renderer/main -i {} -o {} --width {} --height {} -t {} -m'.format(sr, sr.replace('.sr', '.png'), width, height, theme))
    

###########

if generate:
    print("generating animation...")

    import cv2

    pngs = [ os.path.join(folder, x) for x in os.listdir(folder) if x.endswith('.png')]
    pngs = sorted(pngs)

    video = None

    for idx, p in enumerate(pngs):
        print('Generating {} frame'.format(idx))
        f = cv2.imread(p)
        f = cv2.resize(f, (resize_width, resize_height))
        if video == None:
            h, w, l = f.shape
            video = cv2.VideoWriter('{}.avi'.format(avi), cv2.VideoWriter_fourcc(*'XVID'), 30, (w, h))

        if idx == 0:
            for i in range(40):
                video.write(f)
        video.write(f)

    video.release()

print("done")
