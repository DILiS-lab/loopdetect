from matplotlib import rcParams

rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "mathtext.fontset": "cm",
    "text.usetex": True,  # requires LaTeX installation
    "axes.linewidth": 1.2,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    "figure.dpi": 300,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size": 6,
    "xtick.minor.size": 3,
    "ytick.major.size": 6,
    "ytick.minor.size": 3,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
})

colors = ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB']
colors_01 = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499', '#DDDDDD']
colors_02 = ['#332288', '#AA4499', '#DDCC77', '#117733', '#88CCEE', '#44AA99', '#999933', '#CC6677', '#882255', '#DDDDDD']
colors_okabe_ito = ['#000000', '#CC79A7', '#0072B2', '#009E73','#E69F00', '#56B4E9', '#F0E442', '#D55E00']
colors_okabe_ito_2 = ['#D55E00', '#0072B2', '#009E73','#E69F00', '#56B4E9', '#F0E442', '#D55E00', '#CC79A7','#000000']
colors_paul_tol_bright = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']