"""Model parameters and equations."""
PARAMS = dict(
    ft=0.5,               # Free transcription-factor fraction.
    bratio=200.0,         # Exogenous/endogenous mRNA production-rate ratio.
    fn=0.9,               # Inactive RNAP fraction, free or bound to noncoding DNA.
    fr=0.3,               # Free ribosome fraction.
    Nn=3e4,               # RNAP molecules per cell.
    Nr=3e5,               # Ribosomes per cell.
    g=6000.0,             # Number of endogenous genes.
    m_ex=500.0,           # Exogenous mRNAs per gene copy.
    tau_ex=900.0,         # Exogenous mRNA lifetime in seconds.
    Ltx=1000.0,           # Transcribed gene length in nucleotides.
    vtx=30.0,             # Transcription elongation speed in nt/s.
    gamma=0.4 / 3600.0,   # Dilution rate per second (0.4 per hour).
    vtl=10.0,             # Translation elongation speed in aa/s.
    Ttl=24.0,             # Exogenous-protein translation time in seconds.
    p_ex=1.2e6,           # Exogenous proteins per gene copy.
    Mtot=4.0e-12,         # Total cellular protein mass in grams.
    a_ex=27e3 * 1.66053906660e-24,  # Mass of one exogenous protein in grams.
)
PARAMS["Ttx"] = PARAMS["Ltx"] / PARAMS["vtx"]  # Transcription time in seconds.

FT_C, FT_LO, FT_HI = PARAMS["ft"], 0.2, 0.8  # Reference, lower, and upper ft values.
tau_ratio_DAmP_to_wt = 0.1  # DAmP/WT mRNA lifetime and expression ratio.
PGK1_expression_ratio_to_TDH3 = 0.5  # PGK1/TDH3 expression ratio.
CLN2_relative_protein_abundance = 0.05  # Destabilized-GFP abundance relative to GFP.
FLUORESCENCE_TO_PHI = 1.63  # Proteome % per 10^4 fluorescence units.

def c_rnap(p):
    return p["fr"] * p["Ttx"] / (p["Nn"] * p["tau_ex"]) * p["m_ex"]


def c_tf(p):
    return p["bratio"] / p["g"] * (1 - p["ft"]) * p["fn"] * p["fr"]


def c_ribo(p):
    return p["gamma"] * p["Ttl"] / p["Nr"] * p["p_ex"]


def _parameters(**param_overrides):
    unknown = set(param_overrides).difference(PARAMS)
    if unknown:
        names = ", ".join(sorted(unknown))
        raise TypeError(f"unknown model parameter override(s): {names}")
    p = dict(PARAMS)
    p.update(param_overrides)
    if "vtx" in param_overrides or "Ltx" in param_overrides:
        p["Ttx"] = p["Ltx"] / p["vtx"]
    return p


def transcription_cost(**param_overrides):
    p = _parameters(**param_overrides)
    return 100 * (c_rnap(p) + c_tf(p))


def translation_cost(**param_overrides):
    p = _parameters(**param_overrides)
    if "vtl" in param_overrides:
        p["Ttl"] = PARAMS["Ttl"] * PARAMS["vtl"] / p["vtl"]
    return 100 * c_ribo(p)


def total_cost(**param_overrides):
    return transcription_cost(**param_overrides) + translation_cost(**param_overrides)


def tf_fraction(**param_overrides):
    p = _parameters(**param_overrides)
    return c_tf(p) / (c_tf(p) + c_rnap(p))


def cost_per_gene_copy(p):
    return c_rnap(p) + c_tf(p) + c_ribo(p)


def cost_per_proteome_fraction(p):
    beta_m_ex = p["m_ex"] / p["tau_ex"]
    beta_p_ex = p["gamma"] * p["p_ex"] / p["m_ex"]
    return p["gamma"] * p["Mtot"] / p["a_ex"] * (
        (c_tf(p) / beta_m_ex + p["Ttx"] / p["Nn"] * p["fr"])
        / (p["tau_ex"] * beta_p_ex)
        + p["Ttl"] / p["Nr"])


def construct_parameters(name, ft, scenario="weak_TF"):
    p = dict(PARAMS)
    p["ft"] = ft
    if name == "DAmP":
        for key in ("tau_ex", "m_ex", "p_ex"):
            p[key] *= tau_ratio_DAmP_to_wt
    elif name == "PGK1":
        p["m_ex"] *= PGK1_expression_ratio_to_TDH3
        p["p_ex"] *= PGK1_expression_ratio_to_TDH3
        if scenario == "weak_TF":
            p["bratio"] *= PGK1_expression_ratio_to_TDH3
        elif scenario != "weak_RNAP":
            raise ValueError("unknown PGK1 scenario: " + scenario)
    elif name not in ("TDH3", "GFP"):
        raise ValueError("unknown construct: " + name)
    return p


def theory_cost_per_copy(name, ft, scenario="weak_TF"):
    return 100 * cost_per_gene_copy(construct_parameters(name, ft, scenario))


def theory_cost_per_phi(name, ft, scenario="weak_TF"):
    return cost_per_proteome_fraction(construct_parameters(name, ft, scenario))
