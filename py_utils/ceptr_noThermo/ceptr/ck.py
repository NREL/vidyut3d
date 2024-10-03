"""CK routines."""

import ceptr.constants as cc
import ceptr.thermo as cth
import ceptr.writer as cw


def ckawt(fstream, mechanism):
    """Write ckawt."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get atomic weight for all elements"))
    cw.writer(fstream, "void CKAWT" + cc.sym + "( amrex::Real *  awt)")
    cw.writer(fstream, "{")
    cw.writer(fstream, "atomicWeight(awt);")
    cw.writer(fstream, "}")


def ckncf(fstream, mechanism, species_info):
    """Write ckncf."""
    n_elements = mechanism.n_elements
    n_species = species_info.n_species

    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns the elemental composition "))
    cw.writer(fstream, cw.comment("of the speciesi (mdim is num of elements)"))
    cw.writer(fstream, "void CKNCF" + cc.sym + "(int * ncf)")
    cw.writer(fstream, "{")
    cw.writer(fstream, f"int kd = {n_elements}; ")
    cw.writer(fstream, cw.comment("Zero ncf"))
    cw.writer(fstream, f"for (int id = 0; id < kd * {n_species}; ++ id) {{")
    cw.writer(fstream, " ncf[id] = 0; ")
    cw.writer(fstream, "}")

    cw.writer(fstream)
    for sp in species_info.nonqssa_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqssa_species[spec_idx]
        cw.writer(fstream, cw.comment(f"{species.name}"))
        for elem, coef in mechanism.species(sp).composition.items():
            cw.writer(
                fstream,
                f"ncf[ {species_info.ordered_idx_map[sp]} * kd +"
                f" {mechanism.element_index(elem)} ] = {int(coef)}; "
                + cw.comment(f"{elem}"),
            )
        cw.writer(fstream)
    cw.writer(fstream, "}")


def cksyme_str(fstream, mechanism, species_info):
    """Write cksyme."""
    n_elements = mechanism.n_elements
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns the vector of strings of element names"))
    cw.writer(
        fstream,
        "void CKSYME_STR" + cc.sym + "(amrex::Vector<std::string>& ename)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, f"ename.resize({n_elements});")
    for elem in mechanism.element_names:
        cw.writer(
            fstream,
            f'ename[{mechanism.element_index(elem)}] = "{elem}";',
        )
    cw.writer(fstream, "}")


def cksyms_str(fstream, mechanism, species_info):
    """Write cksyms."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns the vector of strings of species names"))
    cw.writer(
        fstream,
        "void CKSYMS_STR" + cc.sym + "(amrex::Vector<std::string>& kname)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, f"kname.resize({n_species});")
    for species in species_info.nonqssa_species_list:
        cw.writer(
            fstream,
            f'kname[{species_info.ordered_idx_map[species]}] = "{species}";',
        )
    cw.writer(fstream, "}")


def ckindx(fstream, mechanism, species_info):
    """Write ckindx."""
    cw.writer(fstream)
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("A few mechanism parameters"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX"
        + cc.sym
        + "(int& mm, int& kk, int& ii, int& nfit)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, f"mm = {mechanism.n_elements};")
    cw.writer(fstream, f"kk = {species_info.n_species};")
    cw.writer(fstream, f"ii = {mechanism.n_reactions};")
    cw.writer(fstream, "nfit = -1; " + cw.comment("Why do you need this anyway ? "))
    cw.writer(fstream, "}")


def ckrp(fstream, mechanism, species_info):
    """Write ckrp."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment(" Returns R, Rc, Patm"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRP"
        + cc.sym
        + "(amrex::Real& ru, amrex::Real& ruc, amrex::Real& pa)",
    )
    cw.writer(fstream, "{")
    cw.writer(
        fstream,
        f" ru  = {(cc.R * cc.ureg.mole * cc.ureg.kelvin / cc.ureg.erg).m:1.14e}; ",
    )
    cw.writer(
        fstream,
        f" ruc = {(cc.Rc * (cc.ureg.mole * cc.ureg.kelvin / cc.ureg.cal)).m:.20f}; ",
    )
    cw.writer(fstream, f" pa  = {cc.Patm:g}; ")
    cw.writer(fstream, "}")

def ckwt(fstream, mechanism, species_info):
    """Write ckwt."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get molecular weight for all species"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT"
        + cc.sym
        + "( amrex::Real wt[])",
    )
    cw.writer(fstream, "{")
    # call molecularWeight
    cw.writer(fstream, "get_mw(wt);")
    cw.writer(fstream, "}")

def ckwc(fstream, mechanism, species_info):
    """Write ckwc."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the production rate for each species"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC"
        + cc.sym
        + "(const amrex::Real T, amrex::Real C[], amrex::Real wdot[], const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)",
    )
    cw.writer(fstream, "{")

    # convert C to SI units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to SI"))
    cw.writer(fstream, f"for (int id = 0; id < {n_species}; ++id) {{")
    cw.writer(fstream, "C[id] *= 1.0e6;")
    cw.writer(fstream, "}")

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "productionRate(wdot, C, T, Te, EN, ener_exch);")

    # convert C and wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, f"for (int id = 0; id < {n_species}; ++id) {{")
    cw.writer(fstream, "C[id] *= 1.0e-6;")
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")

    cw.writer(fstream, "}")


def ckwyp(fstream, mechanism, species_info):
    """Write ckwyp."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns the molar production rate of species"))
    cw.writer(fstream, cw.comment("Given P, T, and mass fractions"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP"
        + cc.sym
        + "(const amrex::Real P, const amrex::Real T,"
        + "const amrex::Real y[], amrex::Real wdot[],"
        + "const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        f"amrex::Real c[{n_species}]; " + cw.comment("temporary storage"),
    )
    cw.writer(fstream, "amrex::Real YOW = 0; ")
    cw.writer(fstream, "amrex::Real PWORT; ")
    cw.writer(fstream)

    # compute inverse of mean molecular weight first (eq 3)
    cw.writer(fstream, cw.comment("Compute inverse of mean molecular wt first"))
    cw.writer(
        fstream,
        f"for (int i = 0; i < {len(species_info.nonqssa_species_list)}; i++)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "YOW += y[i]*imw(i);")
    cw.writer(fstream, "}")

    cw.writer(fstream, cw.comment("PW/RT (see Eq. 7)"))
    cw.writer(
        fstream,
        "PWORT = P/(YOW *"
        f" {(cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m:1.14e} *"
        " T); ",
    )

    cw.writer(fstream, cw.comment("multiply by 1e6 so c goes to SI"))
    cw.writer(fstream, "PWORT *= 1e6; ")

    # now compute conversion
    cw.writer(fstream, cw.comment("Now compute conversion (and go to SI)"))
    cw.writer(
        fstream,
        f"for (int i = 0; i < {len(species_info.nonqssa_species_list)}; i++)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "c[i] = PWORT * y[i]*imw(i);")
    cw.writer(fstream, "}")

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "productionRate(wdot, c, T, Te, EN, ener_exch);")

    # convert wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, f"for (int id = 0; id < {n_species}; ++id) {{")
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckwxp(fstream, mechanism, species_info):
    """Write ckwxp."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns the molar production rate of species"))
    cw.writer(fstream, cw.comment("Given P, T, and mole fractions"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP"
        + cc.sym
        + "(const amrex::Real P, const amrex::Real T,"
        + "const amrex::Real x[], amrex::Real wdot[],"
        + "const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        f"amrex::Real c[{n_species}]; " + cw.comment("temporary storage"),
    )

    cw.writer(
        fstream,
        "amrex::Real PORT = 1e6 *"
        f" P/({(cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m:1.14e} *"
        " T); "
        + cw.comment("1e6 * P/RT so c goes to SI units"),
    )

    # now compute conversion
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute conversion, see Eq 10"))
    cw.writer(fstream, f"for (int id = 0; id < {n_species}; ++id) {{")
    cw.writer(fstream, "c[id] = x[id]*PORT;")
    cw.writer(fstream, "}")

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "productionRate(wdot, c, T, Te, EN, ener_exch);")

    # convert wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, f"for (int id = 0; id < {n_species}; ++id) {{")
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckwyr(fstream, mechanism, species_info):
    """Write ckwyr."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns the molar production rate of species"))
    cw.writer(fstream, cw.comment("Given rho, T, and mass fractions"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR"
        + cc.sym
        + "(const amrex::Real rho, const amrex::Real T, const amrex::Real"
        " y[], "
        " amrex::Real wdot[],"
        + "const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        f"amrex::Real c[{n_species}]; " + cw.comment("temporary storage"),
    )
    cw.writer(fstream)

    # now compute conversion
    cw.writer(fstream, cw.comment("See Eq 8 with an extra 1e6 so c goes to SI"))
    cw.writer(
        fstream,
        f"for (int i = 0; i < {len(species_info.nonqssa_species_list)}; i++)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "c[i] = 1e6 * rho * y[i]*imw(i);")
    cw.writer(fstream, "}")

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("call productionRate"))
    cw.writer(fstream, "productionRate(wdot, c, T, Te, EN, ener_exch);")

    # convert wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, f"for (int id = 0; id < {n_species}; ++id) {{")
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckwxr(fstream, mechanism, species_info):
    """Write ckwxr."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns the molar production rate of species"))
    cw.writer(fstream, cw.comment("Given rho, T, and mole fractions"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR"
        + cc.sym
        + "(const amrex::Real rho, const amrex::Real T, const amrex::Real"
        " x[], "
        " amrex::Real wdot[],"
        + "const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        f"amrex::Real c[{n_species}]; " + cw.comment("temporary storage"),
    )

    cw.writer(
        fstream,
        "amrex::Real XW = 0; " + cw.comment("See Eq 4, 11 in CK Manual"),
    )
    cw.writer(fstream, "amrex::Real ROW; ")

    # compute mean molecular weight first (eq 3)
    cw.writer(fstream, cw.comment("Compute mean molecular wt first"))
    for sp in range(n_species):
        species = species_info.nonqssa_species[sp]
        cw.writer(
            fstream,
            f"XW += x[{species.idx}]*{species.weight:f}; "
            + cw.comment(f"{species.name}"),
        )

    # now compute conversion
    cw.writer(fstream, cw.comment("Extra 1e6 factor to take c to SI"))
    cw.writer(fstream, "ROW = 1e6*rho / XW;")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute conversion, see Eq 11"))
    cw.writer(fstream, f"for (int id = 0; id < {n_species}; ++id) {{")
    cw.writer(fstream, "c[id] = x[id]*ROW;")
    cw.writer(fstream, "}")

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "productionRate(wdot, c, T, Te, EN, ener_exch);")

    # convert wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, f"for (int id = 0; id < {n_species}; ++id) {{")
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckchrg(fstream, self):
    """Write the species unit charge number."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment(" species unit charge number "))
    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void")
    cw.writer(fstream, "CKCHRG(int kcharge[])")
    cw.writer(fstream, "{")
    for i in range(0, self.species_info.n_species):
        species = self.species_info.nonqssa_species[i]
        text = f"kcharge[{i}] = {(int(species.charge)):d};"
        cw.writer(fstream, text + cw.comment(f"{species.name}"))
    cw.writer(fstream, "}")


def ckchrgmass(fstream, species_info):
    """Write the species charge per unit mass."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(fstream, cw.comment(" species charge per unit mass "))
    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void")
    cw.writer(fstream, "CKCHRGMASS(amrex::Real zk[])")
    cw.writer(fstream, "{")
    cw.writer(fstream)
    cw.writer(fstream, f"int kchrg[{n_species}];")
    cw.writer(fstream, "CKCHRG(kchrg);")
    cw.writer(fstream)
    cw.writer(fstream, f"for (int id = 0; id < {n_species}; ++id) {{")
    cw.writer(
        fstream,
        f"zk[id] = {cc.Na:.8e} * {cc.qc:.8e} * kchrg[id] * imw(id);",
    )
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")

# NEED TO DEAL WITH THIS WHEN QSS
def ckinu(fstream, mechanism, species_info, reaction_info, write_sk=False):
    """Write ckinu/skinu."""
    n_reactions = mechanism.n_reactions
    n_gas_reactions = reaction_info.n_reactions
    phase = "surface" if write_sk else "gas"
    function_prefix = "S" if write_sk else "C"
    function_args = (
        "int* /*ki*/, int* /*nu*/" if n_reactions == 0 else "int ki[], int nu[]"
    )

    maxsp = 0

    ns = [0 for _ in range(n_reactions)]
    ki = [[] for _ in range(n_reactions)]
    nu = [[] for _ in range(n_reactions)]

    for orig_idx, _ in reaction_info.idxmap.items():

        # ignore heterogeneous reactions for CKINU and homogeneous reactions for SKINU
        if (phase == "gas" and orig_idx >= n_gas_reactions) or (
            phase == "surface" and orig_idx < n_gas_reactions
        ):
            continue
        # ensure orig_idx is in the range 0, NUM_SURFACE_REACTIONS for SKINU
        if phase == "surface":
            orig_idx -= n_gas_reactions

        reaction = mechanism.reaction(orig_idx)

        for symbol, coefficient in reaction.reactants.items():
            ki[orig_idx].append(species_info.ordered_idx_map[symbol])
            nu[orig_idx].append(-int(coefficient))
        for symbol, coefficient in reaction.products.items():
            ki[orig_idx].append(species_info.ordered_idx_map[symbol])
            nu[orig_idx].append(int(coefficient))

        maxsp = max(maxsp, len(ki[orig_idx]))

    for orig_idx, _ in reaction_info.idxmap.items():
        # ignore heterogeneous reactions for CKINU and homogeneous reactions for SKINU
        if (phase == "gas" and orig_idx >= n_gas_reactions) or (
            phase == "surface" and orig_idx < n_gas_reactions
        ):
            continue
        # ensure orig_idx is in the range 0, NUM_SURFACE_REACTIONS for SKINU
        if phase == "surface":
            orig_idx -= n_gas_reactions

        reaction = mechanism.reaction(orig_idx)

        ns[orig_idx] = len(ki[orig_idx])
        for _ in range(ns[orig_idx], maxsp):
            ki[orig_idx].append(0)
            nu[orig_idx].append(0)

    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            f"Returns a count of {phase} species in a {phase} "
            "reaction, and their indices"
        ),
    )
    cw.writer(fstream, cw.comment("and stoichiometric coefficients. (Eq 50)"))
    cw.writer(
        fstream,
        f"void {function_prefix}KINU"
        + cc.sym
        + f"(const int i, int& nspec, {function_args})",
    )
    cw.writer(fstream, "{")

    if n_reactions > 0:
        str_ns = ",".join(str(x) for x in ns)
        cw.writer(
            fstream,
            f"const int ns[NUM_{phase.upper()}_REACTIONS] =\n     {{{str_ns:s}}};",
        )

        str_ki = ",".join(",".join(str(x) for x in ki[j]) for j in range(n_reactions))
        cw.writer(
            fstream,
            f"const int kiv[NUM_{phase.upper()}_REACTIONS*{maxsp}] =\n    "
            f" {{{str_ki:s}}};",
        )

        str_nu = ",".join(",".join(str(x) for x in nu[j]) for j in range(n_reactions))
        cw.writer(
            fstream,
            f"const int nuv[NUM_{phase.upper()}_REACTIONS*{maxsp}] =\n    "
            f" {{{str_nu:s}}};",
        )

    cw.writer(fstream, "if (i < 1) {")

    cw.writer(fstream, cw.comment("Return max num species per reaction"))
    cw.writer(fstream, f"nspec = {maxsp};")
    cw.writer(fstream, "} else {")

    if n_reactions == 0:
        cw.writer(fstream, "nspec = -1;")
    else:
        cw.writer(fstream, f"if (i > NUM_{phase.upper()}_REACTIONS) {{")
        cw.writer(fstream, "nspec = -1;")
        cw.writer(fstream, "} else {")

        cw.writer(fstream, "nspec = ns[i-1];")
        cw.writer(fstream, "for (int j=0; j<nspec; ++j) {")
        cw.writer(fstream, f"ki[j] = kiv[(i-1)*{maxsp} + j] + 1;")
        cw.writer(fstream, f"nu[j] = nuv[(i-1)*{maxsp} + j];")
        cw.writer(fstream, "}")
        cw.writer(fstream, "}")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")
