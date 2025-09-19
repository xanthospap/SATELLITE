import numpy as np
import pyneb as pn
import satellite.roman as sr


def sum_element_abundancies(element, ionic_abundancies):
    atomic_sums = {}
    for entry in ionic_abundancies:
        if entry["element"] == element:
            spectrum = entry["spectrum"]
            atomic = entry["atomic"]
            if spectrum in atomic_sums:
                t0 = (atomic_sums[spectrum][0] + entry["abundance"]) / 2e0
                t1 = (
                    np.sqrt(
                        (
                            atomic_sums[spectrum][1] * atomic_sums[spectrum][1]
                            + entry["abundance_error"] * entry["abundance_error"]
                        )
                    )
                    / 2e0
                )
                atomic_sums[spectrum] = (t0, t1)
            else:
                atomic_sums[spectrum] = (entry["abundance"], entry["abundance_error"])
    return (
        sum([x[0] for x in atomic_sums.values()]),
        np.sqrt(sum([x[1] * x[1] for x in atomic_sums.values()])) / len(atomic_sums),
        {
            j[0]: j[1]
            for j in zip(
                ["{:}{:}".format(element, sr.roman2int(x)) for x in atomic_sums],
                [(x[0], x[1]) for x in atomic_sums.values()],
            )
        },
    )


def computeAbundancies(fitsd, ionic_abundancies, logger):
    atomic_sums = {}
    for element in list(set([j["element"] for j in fitsd])):
        asum, asum_err, partial_list = sum_element_abundancies(
            element, ionic_abundancies
        )
        atomic_sums.update(partial_list)
    return atomic_sums
