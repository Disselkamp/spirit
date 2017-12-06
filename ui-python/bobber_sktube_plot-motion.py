import pkg_resources
pkg_resources.require("matplotlib==1.5.3")
import numpy as np
import matplotlib.pyplot as plt

# plt.rc('text', usetex=True)
# plt.rc('font', **{'family': 'serif', 'serif': ['charter']})
params = {'text.usetex': True, 'text.latex.preamble': [r'\usepackage{charter}', r'\usepackage{amssymb, amsmath}', r'\usepackage{bm}', r'\usepackage{arevmath}']}
plt.rcParams.update(params)

colors = [(55,126,184),(228,26,28),(255,127,0),(77,175,74),(152,78,163),(255,255,51),(0,0,0)]
#colors = [(228,26,28),(255,127,0),(77,175,74),(152,78,163),(55,126,184),(0,0,0)]
colors = [tuple([x/255. for x in y]) for y in colors]

c = 0
# %\definecolor{color1}{RGB}{55,126,184}   % Blue
# %\definecolor{color2}{RGB}{228,26,28}    % Red
# %\definecolor{color3}{RGB}{255,127,0}    % Orange
# %\definecolor{color4}{RGB}{77,175,74}    % Green
# %\definecolor{color5}{RGB}{152,78,163}   % Purple
# %\definecolor{color6}{RGB}{255,255,51}   % Yellow
# %\definecolor{color01}{RGB}{255,255,255} % White
# %\definecolor{color02}{RGB}{0,0,0}       % Black
# %\definecolor{color03}{RGB}{125,125,125} % Gray


stt_magnitudes = [0.050, 0.060, 0.065, 0.070, 0.075, 0.080, 0.090, 0.100]

for beta in [0.00, 0.02]:
    disy = []
    av = []

    for stt_magnitude in stt_magnitudes:
        pos = np.genfromtxt("../output/12x12x8/beta_%.2f/stt_%.3f/positions_STTmagn%.3f"%(beta, stt_magnitude, stt_magnitude), skip_header=1,)[:, -1]

        N = 12  # system: NxNxz

        dis = np.array([i % N for i in pos])-(pos[0] % N)
        i = 0
        while i < len(pos):
            if pos[i] % N  == 0:
                dis[(i):] += N
                i += 3
            i += 1

        fig1 = plt.figure("beta = %.2f (parallel to x)" % (beta))
        plt.plot(range(0, len(pos)), dis, label="%.3f"%stt_magnitude)
        plt.xlabel("steps (1000 iterations per step)", fontsize=20)
        plt.ylabel("distance [$a$]", fontsize=20)
        plt.legend(loc="best")
        fig1.savefig("position%.2f.eps"%beta, format="eps")

        trav = dis[-100:]

        print beta, stt_magnitude, trav[-1], trav[0]

        N = len(trav)
        time_step = 0.006084
        v = (trav[-1]-trav[0])/(N*1000*time_step)
        av += [v]

        vy = 0
        i = 0
        for i in range(0, len(pos)-1):
            if pos[i] % N != 0 and abs(pos[i+1]-pos[i]) > N/2.:
                vy += 1

        disy += [vy]
    
    print "vt: ", (beta/0.02)
    print "v: ", (av[-2]-av[0])/(stt_magnitudes[-2]-stt_magnitudes[0])

    fig2 = plt.figure(2)

    if beta == 888.00:
        plt.plot(stt_magnitudes, av, "-o", label="perpendicular current", color=colors[c])
    elif beta != 0.16:
        plt.plot(stt_magnitudes, av, "-o", label="$\\beta = %.2f$" % beta, color=colors[c])
        # plt.plot(stt_magnitudes, stt_magnitudes, "--")
    elif beta == 0.16:
        plt.plot(stt_magnitudes[:-1], av[:-1], "-o", label="$\\beta = %.2f$" % beta, color=colors[c])
    c+=1

plt.xlabel(r"$u\mu_s/\gamma$", fontsize=16)
plt.ylabel(r"$\langle v \rangle \mu_s/\gamma $", fontsize=16)
plt.legend(loc="best")
fig2.savefig('skyrmion_motion.eps', format='eps')
