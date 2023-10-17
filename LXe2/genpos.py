import sys

list_scint = []
scint_x_l = 4
scint_y_w = 96
scint_z_h = 1
gap_x = 0.8
gap_z = 1.2
nscintx = 19
nscintz = 19
idx = 0

# colorface = (0., 0.5, 0.5, 0.9)
colorface = (1., 1., 1., 1.)
colorface1 = (1., 0.98, 0., 0.)
colorface2 = (1., 0., 0., 0.)

colorfaces = []
colorfaces1 = []
colorfaces2 = []
for i in range(6):
    colorfaces.append(colorface)
    colorfaces1.append(colorface1)
    colorfaces2.append(colorface2)

coloredges = []
# coloredge1 = (0.4, 0.4, 0.4, 0.9)
coloredge1 = (0., 0.5, 0.5, 0.8)
coloredge2 = (0., 0., 0., 0.)
scintpos = []
zpos = []
for i in range(12):
    if (i>7):
        coloredges.append(coloredge1)
    else:
        coloredges.append(coloredge2)

coloredges2 = []
for i in range(12):
    if (i==1 or i == 3 or i == 5 or i == 7):
        coloredges2.append(coloredge1)
    else:
        coloredges2.append(coloredge2)

for dim in range(2):
    dgapz = -(nscintz+1)*(scint_z_h+gap_z)/2
    for z in range(nscintz):
        if (dim==0):
            if (z%2==0):
                dgapx = -nscintx*(scint_x_l+gap_x)/2
                nscintx_add = nscintx
            else:
                dgapx = -nscintx*(scint_x_l+gap_x)/2-scint_x_l/2
                nscintx_add = nscintx+1
        else:
            if (z%2!=0):
                dgapx = -nscintx*(scint_x_l+gap_x)/2
                nscintx_add = nscintx
            else:
                dgapx = -nscintx*(scint_x_l+gap_x)/2-scint_x_l/2
                nscintx_add = nscintx+1
        for i in range(nscintx_add):
            if (dim==0):
                scintpos.append([dgapx,-scint_y_w/2,dgapz,dim])
                if (i==0):
                    zpos.append(dgapz)
            else:
                scintpos.append([-scint_y_w/2,dgapx,dgapz+(scint_z_h+gap_z)/2,dim])
            dgapx+=(scint_x_l+gap_x)
        dgapz+=(scint_z_h+gap_z)

## Convert to geant4 cordination
for idx,i in enumerate(scintpos):
    if (i[3]==0):
        # scint_x_l,scint_y_w,scint_z_h
        scintpos[idx][0] = scintpos[idx][0]+scint_x_l/2
        scintpos[idx][1] = scintpos[idx][1]+scint_y_w/2
        scintpos[idx][2] = scintpos[idx][2]+scint_z_h/2
    else:
        # scint_y_w,scint_x_l,scint_z_h
        scintpos[idx][0] = scintpos[idx][0]+scint_y_w/2
        scintpos[idx][1] = scintpos[idx][1]+scint_x_l/2
        scintpos[idx][2] = scintpos[idx][2]+scint_z_h/2      
        
outtxt = ""
for idx,i in enumerate(scintpos):
    # if (idx>369):
    #     break
    outtxt+=str(idx)+"\t"
    for j in range(len(i)-1):
        outtxt+="%.3f\t" % i[j]
    if (i[j+1]==0):
        outtxt+= "0\t"+str(i[j+1])+"\n"
    else:
        outtxt+= "90\t"+str(i[j+1])+"\n"
print(outtxt)

outtxtref = ""
k = 0
for idx,i in enumerate(zpos):
    outtxtref+="%d\t%.3f\n" % (k,i-0.01)
    k+=1
with open(sys.argv[1],"w") as f:
    f.write(outtxtref)