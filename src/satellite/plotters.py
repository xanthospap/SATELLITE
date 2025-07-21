import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

defaultPlotOptions = {
    "style_sheet": "default",
    "data_points_color": "blue",
    "data_points_line_width": 1.2,
    "line_color": "blue",
    "line_style": "--",
    # "error_bar_color": "#ff9999",
    "error_bar_color": (1.0, 0.3, 0.3, 0.3),
    "error_bar_width": 2,
    "error_bar_capsize": 3,
    "shaded_error_bars": True,
    "error_bar_alpha": 0.2,
}

# Set y-axis to scientific notation
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-0, 0))  # Control when scientific notation kicks in


def loadPlotOptions(fn=None):
    if not fn:
        return defaultPlotOptions


def parseTotalAbundancies(fn):
    missing_entry = (np.nan, np.nan)
    dct = {}
    omegas = {}
    Us = {}

    with open(fn, "r") as fin:
        for line in fin.readlines():
            if line.lstrip().lower().startswith("slit"):
                l = line.replace("-", "").strip().split()
                columns = [x for x in l[1:]]
            elif line.startswith("Omega"):
                l = line.split()
                for j, v in enumerate(l[1:]):
                    omegas[columns[j]] = v
            elif line.startswith("U "):
                l = line.split()
                for j, v in enumerate(l[1:]):
                    Us[columns[j]] = v
            elif len(line) > 5:
                line = line.replace("/", " ")
                element = line[0:5].strip()
                if element not in dct:
                    dct[element] = {}
                slits = line[5:].rstrip().rstrip(",").split(",")
                for j, slit in enumerate(slits):
                    l = slit.split()
                    no_icf = (float(l[0]), float(l[1]))
                    if "TA" not in dct[element]:
                        dct[element]["TA"] = [no_icf]
                    else:
                        dct[element]["TA"].append(no_icf)
                    icfi = 2
                    while icfi < len(l):
                        icf, val = l[icfi].split(":")
                        if icf not in dct[element]:
                            dct[element][icf] = [(float(val), float(l[icfi + 1]))]
                        else:
                            isz = len(dct[element][icf])
                            if isz != j:
                                while len(dct[element][icf]) == j:
                                    dct[element][icf].append(missing_entry)
                            dct[element][icf].append((float(val), float(l[icfi + 1])))
                        icfi += 2
        else:
            pass
    return columns, dct, omegas, Us


def parseTeNeDat(fn):
    columns = []
    data_te = {}
    data_ne = {}
    with open(fn, "r") as fin:
        for line in fin.readlines():
            if line.lstrip().lower().startswith("slit nr."):
                l = line.replace("-", "").strip().split()
                columns = [x for x in l[2:]]
            elif line.lstrip().lower().startswith("temperature"):
                pass
            else:
                ratio = line[0:45].strip()
                l = line[45:].replace("/", " ").split()
                data_te[ratio] = (
                    [float(l[x]) for x in range(0, len(l), 4)],
                    [float(l[x]) for x in range(1, len(l), 4)],
                )
                data_ne[ratio] = (
                    [float(l[x]) for x in range(2, len(l), 4)],
                    [float(l[x]) for x in range(3, len(l), 4)],
                )
    return columns, data_te, data_ne


def genericParser(fn):
    columns = []
    data = {}
    with open(fn, "r") as fin:
        for line in fin.readlines():
            if line.lstrip().lower().startswith("slit nr."):
                l = line.replace("-", "").strip().split()
                columns = [x for x in l[2:]]
            else:
                l = line.split()
                data[l[0]] = (
                    [float(l[x]) for x in range(1, len(l), 2)],
                    [float(l[x]) for x in range(2, len(l), 2)],
                )
    return columns, data


def genericPlotter(
    columns, data, fnout, y_axis_label, title_keyword, barplot=False, optionsFn=None
):
    plotOptions = loadPlotOptions(optionsFn)

    if plotOptions["style_sheet"] is not None:
        plt.style.use(plotOptions["style_sheet"])

    with PdfPages(fnout) as pdf:
        for k, v in data.items():
            fig, ax = plt.subplots()
            # ax.set_title(f"{title_keyword} {k}")
            x = columns
            y = v[0]
            ey = v[1]
            if not barplot:
                # 1. Scatter plot
                ax.scatter(
                    x,
                    y,
                    facecolors=plotOptions["data_points_color"],
                    edgecolors=plotOptions["error_bar_color"],
                    linewidth=plotOptions["data_points_line_width"],
                    zorder=3,
                )
                # 2. Line plot connecting the points
                ax.plot(
                    x,
                    y,
                    color=plotOptions["line_color"],
                    linestyle=plotOptions["line_style"],
                    zorder=1,
                )
                # 3. Error bars
                if not plotOptions["shaded_error_bars"]:
                    ax.errorbar(
                        x,
                        y,
                        yerr=ey,
                        fmt="o",
                        color="blue",
                        ecolor=plotOptions["error_bar_color"],
                        elinewidth=plotOptions["error_bar_width"],
                        capsize=plotOptions["error_bar_capsize"],
                        zorder=2,
                    )
                else:
                    y = np.array(y)
                    ey = np.array(ey)
                    ax.fill_between(
                        x,
                        y - ey,
                        y + ey,
                        color=plotOptions["error_bar_color"],
                        alpha=plotOptions["error_bar_alpha"],
                        zorder=2,
                    )
            else:
                plt.bar(
                    x,
                    y,
                    yerr=ey,
                    color=plotOptions["data_points_color"],
                    capsize=plotOptions["error_bar_capsize"],
                    ecolor=plotOptions["error_bar_color"],
                    edgecolor="black",
                    zorder=3,
                )

            ax.set_xlabel("Slit Nr.")
            # ax.set_ylabel(f"{y_axis_label}")
            ax.set_ylabel(f"{y_axis_label} {k}")
            ax.grid(True)
            ax.yaxis.set_major_formatter(formatter)

            pdf.savefig(fig)
            plt.close(fig)


def plotTotalAbundancies(fn, fnout=None, optionsFn=None):
    columns, data, _, _ = parseTotalAbundancies(fn)
    plotOptions = loadPlotOptions(optionsFn)

    if plotOptions["style_sheet"] is not None:
        plt.style.use(plotOptions["style_sheet"])

    x = columns
    with PdfPages(fnout) as pdf:
        for k, v in data.items():
            fig, ax = plt.subplots()
            # ax.set_title(f"Total Abundancy for {k}")
            for law, tpls in v.items():
                y = [z[0] for z in tpls]
                ey = [z[1] for z in tpls]
                # 1. Scatter plot
                ax.scatter(
                    x,
                    y,
                    label=law,
                    zorder=3,
                )
                # 2. Line plot connecting the points
                ax.plot(
                    x,
                    y,
                    linestyle=plotOptions["line_style"],
                    zorder=1,
                    label="_nolegend_",
                )
                # 3. Error bars
                if not plotOptions["shaded_error_bars"]:
                    ax.errorbar(
                        x,
                        y,
                        yerr=ey,
                        fmt="o",
                        ecolor=plotOptions["error_bar_color"],
                        elinewidth=plotOptions["error_bar_width"],
                        capsize=plotOptions["error_bar_capsize"],
                        zorder=2,
                        label="_nolegend_",
                    )
                else:
                    y = np.array(y)
                    ey = np.array(ey)
                    ax.fill_between(
                        x,
                        y - ey,
                        y + ey,
                        alpha=plotOptions["error_bar_alpha"],
                        zorder=2,
                        label="_nolegend_",
                    )

            ax.set_xlabel("Slit Nr.")
            ax.set_ylabel(f"Total Abundancy for {k}")
            ax.grid(True)
            ax.yaxis.set_major_formatter(formatter)
            ax.legend()

            pdf.savefig(fig)
            plt.close(fig)


def plotLineIntensities(fn, fnout, barplot=False, optionsFn=None):
    return genericPlotter(*genericParser(fn), fnout, "Intensity of", "", barplot)


def plotLineAbundancies(fn, fnout, barplot=False, optionsFn=None):
    return genericPlotter(*genericParser(fn), fnout, "Abundance of", "", barplot)


def plotLineRatios(fn, fnout, barplot=False, optionsFn=None):
    return genericPlotter(*genericParser(fn), fnout, "Line Ratio", "", barplot)


def PlotTeNeDiagnostics(fn, fnout_te, fnout_ne, barplot=False, optionsFn=None):
    columns, data_te, data_ne = parseTeNeDat(fn)
    genericPlotter(columns, data_te, fnout_te, "Temperature Ratio", "", barplot)
    genericPlotter(columns, data_ne, fnout_ne, "Density Ratio", "", barplot)
