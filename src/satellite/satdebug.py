import pyneb as pn
import re


def debug_obs_labels(obs, max_show=50):
    print("type(obs):", type(obs))
    # these exist on pn.Observation
    print("n_lines:", getattr(obs, "n_lines", None))
    print(
        "n_valid_lines:", getattr(obs, "n_valid_lines", None)
    )  # key signal :contentReference[oaicite:2]{index=2}
    print("unique atoms:", obs.getUniqueAtoms())

    labels = [ln.label for ln in obs.lines]
    print("\nFirst labels:")
    for lab in labels[:max_show]:
        print("  ", lab)

    bad = []
    for lab in labels:
        lab0 = lab.rstrip(
            "e"
        )  # trailing 'e' denotes error column :contentReference[oaicite:3]{index=3}
        if "_" not in lab0:
            bad.append((lab, "missing '_'"))
            continue
        ion, wl = lab0.split("_", 1)
        if ion not in pn.LINE_LABEL_LIST:
            bad.append((lab, f"unknown ion {ion}"))
            continue
        if wl not in pn.LINE_LABEL_LIST[ion]:
            bad.append((lab, f"wl '{wl}' not in LINE_LABEL_LIST[{ion}]"))
    print(f"\nUnrecognized labels: {len(bad)}")
    for item in bad[:max_show]:
        print("  ", item)
