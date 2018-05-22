#!/usr/bin/python
import os

input = 'map.sr'
folder = 'output'
output = 'result'
png = 'result.png'
avi = 'result_3'
steps = '1000'
dt = '0.03'
nu = '1'
rho = '1000'
max_iter = '1000'
tol = '1e-12'

width = '7680'
height = '4320'

os.system("rm ./{}/*".format(folder))

############


print("simulating ...")

output_path = os.path.join(folder, output)
os.system("./main {} {} {} {} {} {} {} {} > test.log".format(input, output_path, steps, dt, nu, rho, max_iter, tol))


###############

print("rendering ...")

mapsr = [ os.path.join(folder, x) for x in os.listdir(folder) if x.endswith('.sr')]

mapsr = sorted(mapsr)
for sr in mapsr:
    os.system('./renderer/main {} {} {} {} 1'.format(sr, sr.replace('.sr', '.png'), width, height))


###########

print("generating animation...")

import cv2

pngs = [ os.path.join(folder, x) for x in os.listdir(folder) if x.endswith('.png')]
pngs = sorted(pngs)

video = None

for p in pngs:
    f = cv2.imread(p)
    f = cv2.resize(f, (3840, 2160))
    if video == None:
        h, w, l = f.shape
        video = cv2.VideoWriter('{}.avi'.format(avi), cv2.VideoWriter_fourcc(*'XVID'), 30, (w, h))

    video.write(f)

video.release()

print("done")
