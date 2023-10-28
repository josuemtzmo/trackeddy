import time

tic = time.time()
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({"font.size": 32})
import importlib
import pdb
import random

import cmocean as cm
import matplotlib.gridspec as gridspec
from pylab import *

import trackeddy
import trackeddy.generator.field_generator as fg
import trackeddy.tracking as ttrack
from trackeddy.geometryfunc import *

importlib.reload(ttrack)
import sys

t = 10
n = 20

xx = linspace(10, 12, 200)
yy = linspace(10, 12, 200)

# print("Generate field")
# gf=fg.Generate_field(0.1,0.1,n,xx,yy,'nrand')
# data = gf.assemble_field(t,'Nint')

x = linspace(10, 12, 300)
y = linspace(10, 12, 300)

data = zeros((t, 300, 300))
for tt in range(t):
    gf = fg.Generate_field(0.05, 0.05, randint(5, n), xx, yy, "int")
    data[tt, :, :] = gf.assemble_field(1)

##

################################################################################
################################################################################
#################################### FLAT ######################################
################################################################################
################################################################################

print("No-Int")
preferences = {"ellipse": 0.85, "eccentricity": 0.85, "gaussian": 0.8}
eddytd = {}
eddytdn = {}

t0 = 0

levels = {"max": data.max(), "min": 0.05, "step": 0.05}
eddytd = trackeddy.tracking.analyseddyzt(
    data,
    x,
    y,
    t0,
    t,
    1,
    levels,
    preferences=preferences,
    areamap="",
    mask="",
    maskopt="forcefit",
    destdir="",
    physics="",
    diagnostics=False,
    plotdata=False,
    pprint=False,
    debug=False,
)

####

levels = {"max": data.min(), "min": -0.05, "step": -0.05}
eddytdn = trackeddy.tracking.analyseddyzt(
    data,
    x,
    y,
    t0,
    t,
    1,
    levels,
    preferences=preferences,
    areamap="",
    mask="",
    maskopt="forcefit",
    destdir="",
    physics="",
    diagnostics=False,
    plotdata=False,
    pprint=False,
    debug=False,
)

pos_f = reconstruct_syntetic(shape(data), x, y, eddytd)
neg_f = reconstruct_syntetic(shape(data), x, y, eddytdn)

f_field = pos_f + neg_f

for tt in range(t0, t):
    f = plt.figure(dpi=300, figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2)
    ax1 = plt.subplot(gs[0])
    ax1.pcolormesh(x, y, data[tt, :, :], vmin=-1, vmax=1, cmap=cm.cm.balance)
    ax1.yaxis.set_major_locator(plt.NullLocator())
    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax2 = plt.subplot(gs[1])
    ax2.pcolormesh(f_field[tt, :, :], vmin=-1, vmax=1, cmap=cm.cm.balance)
    # ax2.contour(f_field[tt,:,:])
    ax2.yaxis.set_major_locator(plt.NullLocator())
    ax2.xaxis.set_major_formatter(plt.NullFormatter())
    # ax1.set_title('Assamble: %03d' % tt)
    plt.show()
    # plt.savefig('plots_n/time_%03d.png' %tt)

m_ke_c = []
m_ke_f = []
m_ke_w = []
m_ke_j = []

for tt in range(shape(data)[0]):
    u_c, v_c = geovelfield(data[tt, :, :], x, y)
    u_f, v_f = geovelfield(f_field[tt, :, :], x, y)
    # u_w,v_w = geovelfield(w_field[tt,:,:],x,y)
    # u_j,v_j = geovelfield(j_field[tt,:,:],x,y)
    ke_c = KE(u_c, v_c)
    ke_f = KE(u_f, v_f)
    # ke_w = KE(u_w,v_w)
    # ke_j = KE(u_j,v_j)
    m_ke_c.append(mean(ke_c))
    m_ke_f.append(mean(ke_f))
    # m_ke_w.append(mean(ke_w))
    # m_ke_j.append(mean(ke_j))

import pandas as pd
import seaborn as sns

figure(dpi=300)
data = np.vstack([m_ke_c, m_ke_f]).T
df = pd.DataFrame(data, columns=[r"$KE_c$", r"$KE_r$"])

sys.exit()

df.to_pickle("./ke_validation_f_n")

sys.exit()

################################################################################
################################################################################
#################################### WAVE ######################################
################################################################################
################################################################################

print("Waves")

amplitude = 1
frequency = 20
phase = 1
waves = zeros(shape(data))

X, Y = meshgrid(x, y)
for t in range(0, t):
    r = X + y / 10
    waves[t, :, :] = 0.3 * sin(r * frequency - t + phase)

wave_data = waves + data

levels = {"max": wave_data.max(), "min": 0.05, "step": 0.05}
eddytd = ttrack.analyseddyzt(
    wave_data,
    x,
    y,
    0,
    t,
    1,
    levels,
    preferences=preferences,
    areamap="",
    mask="",
    maskopt="forcefit",
    destdir="",
    physics="",
    diagnostics=False,
    plotdata=False,
    pprint=False,
)

levels = {"max": wave_data.min(), "min": -0.05, "step": -0.05}
eddytdn = ttrack.analyseddyzt(
    wave_data,
    x,
    y,
    0,
    t,
    1,
    levels,
    preferences=preferences,
    areamap="",
    mask="",
    maskopt="forcefit",
    destdir="",
    physics="",
    diagnostics=False,
    plotdata=False,
    pprint=False,
)

pos_w = reconstruct_syntetic(shape(wave_data), x, y, eddytd)
neg_w = reconstruct_syntetic(shape(wave_data), x, y, eddytdn)

w_field = pos_w + neg_w

for tt in range(t0, t):
    f = plt.figure(dpi=300, figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2)
    ax1 = plt.subplot(gs[0])
    ax1.pcolormesh(x, y, wave_data[tt, :, :], vmin=-1, vmax=1, cmap=cm.cm.balance)
    ax1.yaxis.set_major_locator(plt.NullLocator())
    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax2 = plt.subplot(gs[1])
    ax2.pcolormesh(w_field[tt, :, :], vmin=-1, vmax=1, cmap=cm.cm.balance)
    # ax2.contour(w_field[tt,:,:])
    ax2.yaxis.set_major_locator(plt.NullLocator())
    ax2.xaxis.set_major_formatter(plt.NullFormatter())
    # ax1.set_title('Assamble: %03d' % tt)
    plt.savefig("plots_n/time_w_%03d.png" % tt)

################################################################################
################################################################################
#################################### JETS ######################################
################################################################################
################################################################################

print("Jets")

k_y = 3
phase = 1
k_x = 2
jets = zeros(shape(data))
for t in range(0, t):
    r = Y
    k_y = random.uniform(2, 3)
    phase = random.uniform(0, 1)
    k_x = random.uniform(1, 2)
    amp = 0.3
    jets[t, :, :] = amp * cos((k_y * (k_y * Y + phase + sin(k_x * X - t))))
jet_data = jets + data

levels = {"max": jet_data.max(), "min": 0.05, "step": 0.05}
eddytd = ttrack.analyseddyzt(
    jet_data,
    x,
    y,
    0,
    t,
    1,
    levels,
    preferences=preferences,
    areamap="",
    mask="",
    maskopt="forcefit",
    destdir="",
    physics="",
    diagnostics=False,
    plotdata=False,
    pprint=False,
)

levels = {"max": jet_data.min(), "min": -0.05, "step": -0.05}
eddytdn = ttrack.analyseddyzt(
    jet_data,
    x,
    y,
    0,
    t,
    1,
    levels,
    preferences=preferences,
    areamap="",
    mask="",
    maskopt="forcefit",
    destdir="",
    physics="",
    diagnostics=False,
    plotdata=False,
    pprint=False,
)

pos_f = reconstruct_syntetic(shape(jet_data), x, y, eddytd)
neg_f = reconstruct_syntetic(shape(jet_data), x, y, eddytdn)

j_field = pos_f + neg_f

for tt in range(t0, t):
    f = plt.figure(dpi=300, figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2)
    ax1 = plt.subplot(gs[0])
    ax1.pcolormesh(x, y, jet_data[tt, :, :], vmin=-1, vmax=1, cmap=cm.cm.balance)
    ax1.yaxis.set_major_locator(plt.NullLocator())
    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax2 = plt.subplot(gs[1])
    ax2.pcolormesh(j_field[tt, :, :], vmin=-1, vmax=1, cmap=cm.cm.balance)
    # ax2.contour(w_field[tt,:,:])
    ax2.yaxis.set_major_locator(plt.NullLocator())
    ax2.xaxis.set_major_formatter(plt.NullFormatter())
    # ax1.set_title('Assamble: %03d' % tt)
    plt.savefig("plots_n/time_j_%03d.png" % tt)

################################################################################
################################################################################
##################################### KE #######################################
################################################################################
################################################################################

m_ke_c = []
m_ke_f = []
m_ke_w = []
m_ke_j = []

for tt in range(shape(data)[0]):
    u_c, v_c = geovelfield(data[tt, :, :], x, y)
    u_f, v_f = geovelfield(f_field[tt, :, :], x, y)
    # u_w,v_w = geovelfield(w_field[tt,:,:],x,y)
    # u_j,v_j = geovelfield(j_field[tt,:,:],x,y)
    ke_c = KE(u_c, v_c)
    ke_f = KE(u_f, v_f)
    # ke_w = KE(u_w,v_w)
    # ke_j = KE(u_j,v_j)
    m_ke_c.append(mean(ke_c))
    m_ke_f.append(mean(ke_f))
    # m_ke_w.append(mean(ke_w))
    # m_ke_j.append(mean(ke_j))

################################################################################
################################################################################
#################################### PLOT ######################################
################################################################################
################################################################################

import pandas as pd
import seaborn as sns

figure(dpi=300)
data = np.vstack([m_ke_c, m_ke_f]).T
df = pd.DataFrame(data, columns=[r"$KE_c$", r"$KE_r$"])

df.to_pickle("./ke_validation_f_n")

g1 = sns.jointplot(
    x=r"$KE_c$",
    y=r"$KE_r$",
    data=df,
    kind="kde",
    cmap="Blues",
    joint_kws={"shade_lowest": False},
    fontsize=32,
)

lims = [100, 0]
g1.ax_joint.plot(lims, lims, "--k")

res = stats.theilslopes(df[r"$KE_r$"].values, df[r"$KE_c$"].values, 0.95)

lnr2 = res[1] + res[2] * range(100)
lnr3 = res[1] + res[3] * range(100)
g1.ax_joint.fill_between(range(100), lnr2, lnr3, facecolor="b", alpha=0.5)

r = res[0]
x0 = 0
y0 = res[1] + res[0] * x0
x1 = 100
y1 = res[1] + res[0] * x1

g1.ax_joint.plot([x0, x1], [y0, y1], "-.b")
g1.ax_joint.text(60, 20, r"R = %.2f" % r, color="b")
g1.ax_marg_x.set_xlim(0, 100)
g1.ax_marg_y.set_ylim(0, 100)
print("estimate flat: ", mean([abs(y0 / 100), abs(1 - y1 / 100)]))
g1.ax_joint.legend_.remove()
plt.savefig("e_vs_e_n.png")

figure(dpi=300)
data = np.vstack([m_ke_c, m_ke_w]).T
df = pd.DataFrame(data, columns=[r"$KE_c$", r"$KE_r$"])
g1 = sns.jointplot(
    x=r"$KE_c$",
    y=r"$KE_r$",
    data=df,
    kind="kde",
    cmap="Blues",
    joint_kws={"shade_lowest": False},
    fontsize=32,
)

lims = [100, 0]
g1.ax_joint.plot(lims, lims, "--k")

res = stats.theilslopes(df[r"$KE_r$"].values, df[r"$KE_c$"].values, 0.95)

lnr2 = res[1] + res[2] * range(100)
lnr3 = res[1] + res[3] * range(100)
g1.ax_joint.fill_between(range(100), lnr2, lnr3, facecolor="b", alpha=0.5)

r = res[0]
x0 = 0
y0 = res[1] + res[0] * x0
x1 = 100
y1 = res[1] + res[0] * x1

g1.ax_joint.plot([x0, x1], [y0, y1], "-.b")
g1.ax_joint.text(60, 20, r"R = %.2f" % r, color="b")
g1.ax_marg_x.set_xlim(0, 100)
g1.ax_marg_y.set_ylim(0, 100)
print("estimate sin: ", mean([abs(y0 / 100), abs(1 - y1 / 100)]))
g1.ax_joint.legend_.remove()
plt.savefig("w_vs_e_n.png")
# df.to_pickle('./ke_validation_w_n')

figure(dpi=300)
data = np.vstack([m_ke_c, m_ke_j]).T
df = pd.DataFrame(data, columns=[r"$KE_c$", r"$KE_r$"])
g1 = sns.jointplot(
    x=r"$KE_c$",
    y=r"$KE_r$",
    data=df,
    kind="kde",
    cmap="Blues",
    joint_kws={"shade_lowest": False},
    fontsize=32,
)

lims = [100, 0]
g1.ax_joint.plot(lims, lims, "--k")

res = stats.theilslopes(df[r"$KE_r$"].values, df[r"$KE_c$"].values, 0.95)

lnr2 = res[1] + res[2] * range(100)
lnr3 = res[1] + res[3] * range(100)
g1.ax_joint.fill_between(range(100), lnr2, lnr3, facecolor="b", alpha=0.5)

r = res[0]
x0 = 0
y0 = res[1] + res[0] * x0
x1 = 100
y1 = res[1] + res[0] * x1

g1.ax_joint.plot([x0, x1], [y0, y1], "-.b")
g1.ax_joint.text(60, 20, r"R = %.2f" % r, color="b")
g1.ax_marg_x.set_xlim(0, 100)
g1.ax_marg_y.set_ylim(0, 100)
print("estimate jet: ", mean([abs(y0 / 100), abs(1 - y1 / 100)]))
g1.ax_joint.legend_.remove()
plt.savefig("j_vs_e_n.png")
# df.to_pickle('./ke_validation_j_n')


# for ii in range(0,30):
#      plt.figure()
#      plt.pcolormesh(af[ii])
#      plt.savefig('%03d.png' %ii)
#      plt.show()

toc = time.time()

print("######## ELAPSED TIME: ###########")
print("######## %2f s ###########" % (toc - tic))
