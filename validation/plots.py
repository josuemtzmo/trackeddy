plt.contourf(ssh)
plt.colorbar()
plt.savefig('out.png')


plt.plot(CONTS[0][0][:,0],CONTS[0][0][:,1])
plt.plot(CONTS[0][1][:,0],CONTS[0][1][:,1])
plt.plot(CONTS[0][2][:,0],CONTS[0][2][:,1])
plt.plot(CONTS[0][3][:,0],CONTS[0][3][:,1])
plt.plot(CONTS[0][4][:,0],CONTS[0][4][:,1])
plt.plot(CONTS[0][5][:,0],CONTS[0][5][:,1])

plt.plot(CONTS[0][0][:,0],CONTS[0][0][:,1])
plt.plot(CONTS[0][1][:,0],CONTS[0][1][:,1])
plt.plot(CONTS[0][2][:,0],CONTS[0][2][:,1])

plt.contourf(lon,lat,ssh,vmin=-1,vmax=1)
plt.colorbar()
plt.contour(lon,lat,ssh,'-k',levels=[0.2,np.inf])
# plt.plot(CONTeach[:,0],CONTeach[:,1],'*-')
plt.plot(ellipse['ellipse'][0],ellipse['ellipse'][1],'r')
plt.savefig('cont.png')
plt.close()

x=np.linspace(10,12,300)
y=np.linspace(10,12,300)
X,Y=np.meshgrid(x,y)
coords1 = [X,Y,eddz['PositionExtreme'][0][2],x[int(eddz['PositionExtreme'][0][3])],y[int(eddz['PositionExtreme'][0][4])]]
gaus1 = twoD_Gaussian(coords1, *eddz['2DGaussianFit'][0][0])
g1=gaus1.reshape((300,300))

coords2=[X,Y,eddz['PositionExtreme'][1][2],x[int(eddz['PositionExtreme'][1][3])],y[int(eddz['PositionExtreme'][1][4])]]
gaus2 = twoD_Gaussian(coords2, *eddz['2DGaussianFit'][0][0])
g2=gaus2.reshape((300,300))

plt.pcolormesh(g1+g2)
plt.savefig('gauss.png')


x=np.linspace(10,12,300)
y=np.linspace(10,12,300)
X,Y=np.meshgrid(x,y)
key='eddyn_0'
counter=0
coords1 = [X,Y,eddytdn[key]['position_maxvalue'][counter][2],eddytdn[key]['position_maxvalue'][counter][0],eddytdn[key]['position_maxvalue'][counter][1]]
gaus1 = twoD_Gaussian(coords1, *eddytdn[key]['2dgaussianfit'][0])
g1=gaus1.reshape((300,300))

key='eddyn_1'
coords2 = [X,Y,eddytdn[key]['position_maxvalue'][counter][2],eddytdn[key]['position_maxvalue'][counter][0],eddytdn[key]['position_maxvalue'][counter][1]]
gaus2 = twoD_Gaussian(coords2, *eddytdn[key]['2dgaussianfit'][0])
g2=gaus2.reshape((300,300))

plt.pcolormesh(g1+g2)
plt.savefig('gauss_t.png')


gaus2 =twoD_Gaussian(maxposition, *curvefit)
g2=gaus2.reshape((300,300))

plt.pcolormesh(g2)
plt.savefig('gauss_r.png')

plt.pcolormesh(fieldfit[0,:,:])
plt.contour(fieldfit[0,:,:])
plt.contour(fittedcurve.reshape(len(lat),len(lon)))
plt.colorbar()
plt.savefig('recons_bug.png')
plt.close()


gaus2 =twoD_Gaussian(fix, *gausssianfitp)
g2=gaus2.reshape((75,75))

plt.pcolormesh(g2)
plt.savefig('gauss_r.png')

print("| ", areastatus['ellipse'] < areastatus['check'] and areastatus['contour'] < areastatus['check'] and  gaussarea[0]," | ", checke ,"| ",eccen <= preferences['eccentricity'] and fiteccen <= preferences['eccentricity'] ," | ", checkM and checkm)

x=np.linspace(10,12,300)
y=np.linspace(10,12,300)
pos=reconstruct_syntetic((3,300,300),x,y,eddydt)

plt.pcolormesh(pos[2,:,:])
plt.colorbar()
plt.savefig('recons.png')

f, (ax1, ax2,ax3) = plt.subplots(1, 3,figsize=(25,7),sharey=True)
ax1.pcolormesh(values[0],values[1],varm,vmin=-1,vmax=1)
ax1.set_title('Original Field')
ax1.axis('equal')
ax2.pcolormesh(values[0], values[1],fittedata,vmin=-1,vmax=1)
ax2.set_title('2D Gauss Fit')
ax2.axis('equal')
ax3.pcolormesh(values[0], values[1],varm-fittedata,vmin=-1,vmax=1,cmap=cm.cm.balance)
ax3.axis('equal')
ax3.set_title('Difference between Fit & Original')
plt.savefig('fit.png')

res = minimize(gaussian2Dresidual, initial_guess,args=(coords,varm),method='SLSQP',options={'xtol': 1e-12, 'disp': False})
fitdict = res.x
fitted_curve = twoD_Gaussian(coords, *fitdict)
fittedata=fitted_curve.reshape(len(values[1]), len(values[0]))


# [[key,item['2dgaussianfit']] for key,item in eddytd.items()]
# [['eddyn_0', array([[ 0.09007907,  0.10030357, -0.63405094,  0.        ,  0.        ,
#          0.        ],
#        [-0.00155164, -0.00155203,  0.        ,  0.        ,  0.        ,
#          0.        ]])], ['eddyn_1', array([[ 0.08936764,  0.10073398, -0.63314404,  0.        ,  0.        ,
#          0.        ]])], ['eddyn_2', array([[-0.00155164, -0.00155203,  0.        ,  0.        ,  0.        ,
#          0.        ]])], ['eddyn_3', array([[ -9.86395841e-02,  -9.97461553e-02,   1.22649918e-09,
#           0.00000000e+00,   0.00000000e+00,   0.00000000e+00]])], ['eddyn_5', array([[ -9.42872770e-02,  -1.04469383e-01,  -4.28372171e-10,
#           0.00000000e+00,   0.00000000e+00,   0.00000000e+00]])]]


#array([-0.00155164, -0.00155203,  0.        ,  0.        ,  0.        ,  0.        ])
#[0.24077341334104713, 0.24077330578356948, 0, 0, 0, 0]

