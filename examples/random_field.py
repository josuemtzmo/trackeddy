import time

tic = time.time()
import matplotlib

matplotlib.use("Agg")
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


t = 1000
n = 13

xx = linspace(10, 12, 200)
yy = linspace(10, 12, 200)

# print("Generate field")
# gf=fg.Generate_field(0.1,0.1,n,xx,yy,'Nint')
# data = gf.assemble_field(t)

data = zeros((t, 300, 300))
for tt in range(t):
    print(tt)
    gf = fg.Generate_field(0.1, 0.1, randint(5, 15), xx, yy, "Nint")
    data[tt, :, :] = gf.assemble_field(1)

##

x = linspace(10, 12, 300)
y = linspace(10, 12, 300)

################################################################################
################################################################################
#################################### FLAT ######################################
################################################################################
################################################################################

preferences = {"ellipse": 0.85, "eccentricity": 0.85, "gaussian": 0.8}
eddytd = {}
eddytdn = {}

t0 = 0
t = 1000

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
    pprint=True,
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
    pprint=True,
    debug=False,
)

pos_f = reconstruct_syntetic(shape(data), x, y, eddytd)
neg_f = reconstruct_syntetic(shape(data), x, y, eddytdn)

f_field = pos_f + neg_f

for tt in range(t0, t):
    f = plt.figure()
    gs = gridspec.GridSpec(2, 1)
    ax1 = plt.subplot(gs[0])
    ax1.pcolormesh(x, y, data[tt, :, :], vmin=-1, vmax=1, cmap=cm.cm.balance)
    ax2 = plt.subplot(gs[1])
    ax2.pcolormesh(f_field[tt, :, :], vmin=-1, vmax=1, cmap=cm.cm.balance)
    ax2.contour(f_field[tt, :, :])
    ax1.set_title("Assamble: %03d" % tt)
    plt.savefig("time_%03d.png" % tt)

################################################################################
################################################################################
#################################### WAVE ######################################
################################################################################
################################################################################

amplitude = 1
frequency = 20
phase = 1
waves = zeros(shape(data))

X, Y = meshgrid(x, y)
for t in range(0, t):
    r = X + y / 10
    waves[t, :, :] = 0.3 * sin(r * frequency - t + phase)

wave_data = waves + data

levels = {"max": data.max(), "min": 0.05, "step": 0.05}
eddytd = ttrack.analyseddyzt(
    data,
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
    pprint=True,
)

levels = {"max": data.min(), "min": -0.05, "step": -0.05}
eddytdn = ttrack.analyseddyzt(
    data,
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
    pprint=True,
)

pos_w = reconstruct_syntetic(shape(data), x, y, eddytd)
neg_w = reconstruct_syntetic(shape(data), x, y, eddytdn)

w_field = pos_w + neg_w

################################################################################
################################################################################
#################################### JETS ######################################
################################################################################
################################################################################

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

levels = {"max": data.max(), "min": 0.05, "step": 0.05}
eddytd = ttrack.analyseddyzt(
    data,
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
    pprint=True,
)

levels = {"max": data.min(), "min": -0.05, "step": -0.05}
eddytdn = ttrack.analyseddyzt(
    data,
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
    pprint=True,
)

pos_f = reconstruct_syntetic(shape(data), x, y, eddytd)
neg_f = reconstruct_syntetic(shape(data), x, y, eddytdn)

j_field = pos_f + neg_f

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
    u_w, v_w = geovelfield(w_field[tt, :, :], x, y)
    u_j, v_j = geovelfield(j_field[tt, :, :], x, y)
    ke_c = KE(u_c, v_c)
    ke_f = KE(u_f, v_f)
    ke_w = KE(u_w, v_w)
    ke_j = KE(u_j, v_j)
    m_ke_c.append(mean(ke_c))
    m_ke_f.append(mean(ke_f))
    m_ke_w.append(mean(ke_w))
    m_ke_j.append(mean(ke_j))

################################################################################
################################################################################
#################################### PLOT ######################################
################################################################################
################################################################################

import pandas as pd
import seaborn as sns
from scipy.stats import linregress, spearmanr

figure(dpi=300)
data = np.vstack([m_ke_c, m_ke_f]).T
df = pd.DataFrame(data, columns=[r"$KE_c$", r"$KE_r$"])
g1 = sns.jointplot(
    x=r"$KE_c$",
    y=r"$KE_r$",
    data=df,
    kind="kde",
    cmap="Blues",
    joint_kws={"shade_lowest": False},
)

lims = [100, 0]
g1.ax_joint.plot(lims, lims, "--k")

s, i, r, p, std = linregress(m_ke_c, m_ke_f)

x0 = 0
y0 = s * x0 + i
x1 = 100
y1 = s * x1 + i

g1.ax_joint.plot([x0, x1], [y0, y1], "-.b")
g1.ax_joint.text(60, 20, r"R = %2f" % r, color="b")
g1.ax_marg_x.set_xlim(0, 100)
g1.ax_marg_y.set_ylim(0, 100)
print("estimate flat: ", mean([abs(y0 / 100), abs(1 - y1 / 100)]))
plt.savefig("e_vs_e.png")

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
)

lims = [100, 0]
g1.ax_joint.plot(lims, lims, "--k")

s, i, r, p, std = linregress(m_ke_c, m_ke_w)

x0 = 0
y0 = s * x0 + i
x1 = 100
y1 = s * x1 + i

g1.ax_joint.plot([x0, x1], [y0, y1], "-.b")
g1.ax_joint.text(60, 20, r"R = %2f" % r, color="b")
g1.ax_marg_x.set_xlim(0, 100)
g1.ax_marg_y.set_ylim(0, 100)
print("estimate sin: ", mean([abs(y0 / 100), abs(1 - y1 / 100)]))
plt.savefig("w_vs_e.png")

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
)

lims = [100, 0]
g1.ax_joint.plot(lims, lims, "--k")

s, i, r, p, std = linregress(m_ke_c, m_ke_j)

x0 = 0
y0 = s * x0 + i
x1 = 100
y1 = s * x1 + i

g1.ax_joint.plot([x0, x1], [y0, y1], "-.b")
g1.ax_joint.text(60, 20, r"R = %2f" % r, color="b")
g1.ax_marg_x.set_xlim(0, 100)
g1.ax_marg_y.set_ylim(0, 100)
print("estimate jet: ", mean([abs(y0 / 100), abs(1 - y1 / 100)]))
plt.savefig("j_vs_e.png")


# for ii in range(0,30):
#      plt.figure()
#      plt.pcolormesh(af[ii])
#      plt.savefig('%03d.png' %ii)
#      plt.show()

toc = time.time()

print("######## ELAPSED TIME: ###########")
print("######## %2f s ###########" % (toc - tic))
