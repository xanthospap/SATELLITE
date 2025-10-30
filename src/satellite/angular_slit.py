from satellite import specific_slit


def angularOptionsToSlit(ang_options, idx):
    slit = {}
    slit["PA"] = ang_options["start_angle"] + idx * ang_options["step_angle"]
    slit["w"] = ang_options["w"]
    slit["h"] = ang_options["h"]
    slit["x"] = ang_options["x"]
    slit["y"] = ang_options["y"] + ang_options["h"] // 2
    return slit


def angular_slit_analysis(
    fitsd: list,
    slit: list,
    ratios: list,
    density_diagnostics: list,
    tempterature_diagnostics: list,
    ext_law: str,
    pn_atomic_data: str,
    monte_carlo_fake_obs: int,
    max_nan_in_diagnostics_allowed_percentage: int,
    pn_rv: float,
    energy_parameter: float,
    intensities_out: str,
    ratios_out: str,
    diagnostics_out: str,
    abundancies_out: str,
    total_abundancies_out: str,
    corners_out: str,
    logger,
):
    # slits should be a one-element list, of type:
    # {"start_angle":, "stop_angle":, "step_angle":, "w":, "h":, "x":, "y":}
    assert len(slit) == 1
    slit = slit[0]
    assert slit["start_angle"] >= 0e0 and slit["start_angle"] <= 360e0
    assert slit["stop_angle"] > slit["start_angle"]
    assert slit["step_angle"] >= 0e0
    num_slits = int(
        float(slit["stop_angle"] - slit["start_angle"]) / slit["step_angle"]
    )

    slits = [angularOptionsToSlit(slit, j) for j in range(num_slits + 1)]

    # just call specific_slit with the constructed slits ...
    # Specific Slit Analysis
    return specific_slit.specific_slit_analysis(
        fitsd,
        slits,
        ratios,
        density_diagnostics,
        tempterature_diagnostics,
        ext_law,
        pn_atomic_data,
        monte_carlo_fake_obs,
        max_nan_in_diagnostics_allowed_percentage,
        pn_rv,
        energy_parameter,
        intensities_out,
        ratios_out,
        diagnostics_out,
        abundancies_out,
        total_abundancies_out,
        corners_out,
        logger,
    )
