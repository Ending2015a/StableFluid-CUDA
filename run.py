#!/usr/bin/python
import os

input = 'map.sr'
folder = 'output'
output = 'result'
png = 'result.png'
steps = '100'
dt = '0.01'
nu = '1'
rho = '1000'
max_iter = '1000'
tol = '1e-12'

width = '1280'
height = '720'

os.system("rm ./{}/*".format(folder))

############


print("simulating ...")

output_path = os.path.join(folder, output)
os.system("./main {} {} {} {} {} {} {} {}".format(input, output_path, steps, dt, nu, rho, max_iter, tol))


###############

print("rendering ...")

mapsr = [ os.path.join(folder, x) for x in os.listdir(folder) if x.endswith('.sr')]

mapsr = sorted(mapsr)
for sr in mapsr:
    os.system('./renderer/main {} {} {} {}'.format(sr, sr.replace('.sr', '.png'), width, height))


###########

print("generating animation...")

import cv2

pngs = [ os.path.join(folder, x) for x in os.listdir(folder) if x.endswith('.png')]
pngs = sorted(pngs)

video = None

for p in pngs:
    f = cv2.imread(p)
    if video == None:
        h, w, l = f.shape
        video = cv2.VideoWriter('{}.avi'.format(output), cv2.VideoWriter_fourcc(*'XVID'), 30, (w, h))

    video.write(f)

video.release()

print("done")
