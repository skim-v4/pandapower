def calc_sc():

    if use_pre_fault_voltage:
        # sets initial voltage magnitude and angle to "results"
        init_vm_pu = init_va_degree = "results"

        # trafo model for SC must match the model for PF calculation
        trafo_model = net._options["trafo_model"]

    # stops type C fault calculations on networks with multiple buses and sgen
    if not isinstance(bus, Number) and len(net.sgen.query("in_service")) > 0:
        raise NotImplementedError("Short-circuit with Type C method and sgen is only implemented for a single bus")

    # uses nominal bus voltages and trafo model "pi"
    else:
        init_vm_pu = init_va_degree = "flat"
        trafo_model = "pi"

def _calc_sc_1ph(net, bus):

    # a generator for power injection and another for extraction is created for each "dcline" in net for DC modeling
    _add_auxiliary_elements(net)

    # pos. seq bus impedance
    ppc_1, ppci_1 = _init_ppc(net)

    # Create k updated ppci_1
    ppci_bus = _get_is_ppci_bus(net, bus)
    _, ppci_1, _ = _create_k_updated_ppci(net, ppci_1, ppci_bus=ppci_bus)
    _calc_ybus(ppci_1)

    # input for negative sequence is same as the positive sequence
    ppc_2 = copy.deepcopy(ppc_1)
    ppci_2 = copy.deepcopy(ppci_1)

    # avoid the calculation of Ybus if the fault isn't type C
    if net._options.get("use_pre_fault_voltage", False):
        _add_load_sc_impedances_ppc(net, ppc_1)  # add SC impedances for sgens and loads
        ppci_1 = _ppc2ppci(ppc_1, net)
        _, ppci_1, _ = _create_k_updated_ppci(net, ppci_1, ppci_bus=ppci_bus)
        _calc_ybus(ppci_1)

        _add_load_sc_impedances_ppc(net, ppc_2, relevant_elements=("load",))  # add SC impedances for loads
        ppci_2 = _ppc2ppci(ppc_2, net)
        _calc_ybus(ppci_2)

    # zero seq bus impedance
    ppc_0, ppci_0 = _pd2ppc_zero(net, ppc_1['branch'][:, K_ST])
    _calc_ybus(ppci_0)

    if net["_options"]["inverse_y"]:
        _calc_zbus(net, ppci_0)
        _calc_zbus(net, ppci_1)
        _calc_zbus(net, ppci_2)
    else:
        # Factorization Ybus once
        ppci_0["internal"]["ybus_fact"] = factorized(ppci_0["internal"]["Ybus"].tocsc())
        ppci_1["internal"]["ybus_fact"] = factorized(ppci_1["internal"]["Ybus"].tocsc())
        ppci_2["internal"]["ybus_fact"] = factorized(ppci_2["internal"]["Ybus"].tocsc())

    _calc_rx(net, ppci_1, ppci_bus)
    _add_kappa_to_ppc(net, ppci_1)  # todo add kappa only to ppci_1?

    _calc_rx(net, ppci_0, ppci_bus)
    _calc_rx(net, ppci_2, ppci_bus)

    _calc_ikss_1ph(net, ppci_0, ppci_1, ppci_2, ppci_bus)
    # from here on, the V_ikss in ppci_0, ppci_1, ppci_2 are in phase frame!

    if net._options["branch_results"]:
        _calc_branch_currents_complex(net, ppci_0, ppci_bus)
        _calc_branch_currents_complex(net, ppci_1, ppci_bus)
        _calc_branch_currents_complex(net, ppci_2, ppci_bus)

    ppc_0 = _copy_results_ppci_to_ppc(ppci_0, ppc_0, "sc")
    ppc_1 = _copy_results_ppci_to_ppc(ppci_1, ppc_1, "sc")
    ppc_2 = _copy_results_ppci_to_ppc(ppci_2, ppc_2, "sc")
    _extract_results(net, ppc_0, ppc_1, ppc_2, bus=bus)
    _clean_up(net)

def _init_ppc(net):

    # checks all columns in trafo and gen are populated, fills missing values with NaN
    _check_sc_data_integrity(net)

    ppc, _ = _pd2ppc(net)

    # Init the required columns to nan
    ppc["bus"][:, [K_G, K_SG, V_G, PS_TRAFO_IX, GS_P, BS_P, GS_GEN, BS_GEN]] = np.nan
    ppc["branch"][:, [K_T, K_ST]] = np.nan

    # Add parameter K into ppc
    _add_kt(net, ppc)
    _add_gen_sc_z_kg_ks(net, ppc)
    _add_sgen_sc_z(net, ppc)
    _add_ward_sc_z(net, ppc)

    ppci = _ppc2ppci(ppc, net)

    return ppc, ppci

def _pd2ppc(net, sequence=None):
    """
    Converter Flow:
        1. Create an empty pypower data structure
        2. Calculate loads and write the bus matrix
        3. Build the gen (Infeeder)- Matrix
        4. Calculate the line parameter and the transformer parameter,
           and fill it in the branch matrix.
           Order: 1st: Line values, 2nd: Trafo values
        5. if opf: make opf objective (gencost)
        6. convert internal ppci format for pypower powerflow /
        opf without out of service elements and rearanged buses

    INPUT:
        **net** - The pandapower format network
        **sequence** - Used for three phase analysis
        ( 0 - Zero Sequence
          1 - Positive Sequence
          2 - Negative Sequence
        )

    OUTPUT:
        **ppc** - The simple matpower format network. Which consists of:
                  ppc = {
                        "baseMVA": 1., *float*
                        "version": 2,  *int*
                        "bus": np.array([], dtype=float),
                        "branch": np.array([], dtype=np.complex128),
                        "gen": np.array([], dtype=float),
                        "gencost" =  np.array([], dtype=float), only for OPF
                        "internal": {
                              "Ybus": np.array([], dtype=np.complex128)
                              , "Yf": np.array([], dtype=np.complex128)
                              , "Yt": np.array([], dtype=np.complex128)
                              , "branch_is": np.array([], dtype=bool)
                              , "gen_is": np.array([], dtype=bool)
                              }
        **ppci** - The "internal" pypower format network for PF calculations

    """
    # select elements in service (time consuming, so we do it once)
    net["_is_elements"] = aux._select_is_elements_numba(net, sequence=sequence)

    # Gets network configurations
    mode = net["_options"]["mode"]
    check_connectivity = net["_options"]["check_connectivity"]
    calculate_voltage_angles = net["_options"]["calculate_voltage_angles"]

    ppc = _init_ppc(net, mode=mode, sequence=sequence)

    # generate ppc['bus'] and the bus lookup
    _build_bus_ppc(net, ppc, sequence=sequence)

    if sequence == 0:
        from pandapower.pd2ppc_zero import _add_ext_grid_sc_impedance_zero, _build_branch_ppc_zero
        # Adds external grid impedance for 3ph and sc calculations in ppc0
        _add_ext_grid_sc_impedance_zero(net, ppc)
        # Calculates ppc0 branch impedances from branch elements
        _build_branch_ppc_zero(net, ppc)
    else:
        # Calculates ppc1/ppc2 branch impedances from branch elements
        _build_branch_ppc(net, ppc)

    _build_tcsc_ppc(net, ppc, mode)
    _build_svc_ppc(net, ppc, mode)
    _build_ssc_ppc(net, ppc, mode)

    # Adds P and Q for loads / sgens in ppc['bus'] (PQ nodes)
    if mode == "sc":
        _add_ext_grid_sc_impedance(net, ppc)
        # Generator impedance are seperately added in sc module
        _add_motor_impedances_ppc(net, ppc)

    else:
        _calc_pq_elements_and_add_on_ppc(net, ppc, sequence=sequence)
        # adds P and Q for shunts, wards and xwards (to PQ nodes)
        _calc_shunts_and_add_on_ppc(net, ppc)

    # adds auxilary buses for open switches at branches
    _switch_branches(net, ppc)

    # Adds auxilary buses for in service lines with out of service buses.
    # Also deactivates lines if they are connected to two out of service buses
    _branches_with_oos_buses(net, ppc)

    if check_connectivity:
        if sequence in [None, 1, 2]:
            # sets islands (multiple isolated nodes) out of service
            if "opf" in mode:
                net["_isolated_buses"], _, _ = aux._check_connectivity_opf(ppc)
            else:
                net["_isolated_buses"], _, _ = aux._check_connectivity(ppc)
            net["_is_elements_final"] = aux._select_is_elements_numba(net,
                                                                      net._isolated_buses, sequence)
        else:
            ppc["bus"][net._isolated_buses, BUS_TYPE] = NONE
        net["_is_elements"] = net["_is_elements_final"]
    else:
        # sets buses out of service, which aren't connected to branches / REF buses
        aux._set_isolated_buses_out_of_service(net, ppc)

    _build_gen_ppc(net, ppc)

    if "pf" in mode:
        _check_for_reference_bus(ppc)

    aux._replace_nans_with_default_limits(net, ppc)

    # generates "internal" ppci format (for powerflow calc)
    # from "external" ppc format and updates the bus lookup
    # Note: Also reorders buses and gens in ppc
    ppci = _ppc2ppci(ppc, net)

    if mode == "pf":
        # check if any generators connected to the same bus have different voltage setpoints
        _check_voltage_setpoints_at_same_bus(ppc)
        if calculate_voltage_angles:
            _check_voltage_angles_at_same_bus(net, ppci)

    if mode == "opf":
        # make opf objective
        ppci = _make_objective(ppci, net)

    return ppc, ppci

def _select_is_elements_numba(net, isolated_nodes=None, sequence=None):
# is missing sgen_controllable and load_controllable
# filters out of service or isolated elements and buses to be ignored for sc calculations

    # retrieves the max index of buses
    max_bus_idx = np.max(net["bus"].index.values)

    # initialize a boolean array with length of max bus index + 1, all values 0 (out of service)
    bus_in_service = np.zeros(max_bus_idx + 1, dtype=bool)

    # sets each bufs =s index to its respective bhs////9
    bus_in_service[net["bus"].index.values] = net["bus"]["in_service"].values.astype(bool)

    # only executed if there are isolated nodes
    if isolated_nodes is not None and len(isolated_nodes) > 0:
        ppc = net["_ppc"] if sequence is None else net["_ppc%s" % sequence]
        ppc_bus_isolated = np.zeros(ppc["bus"].shape[0], dtype=bool)
        ppc_bus_isolated[isolated_nodes] = True
        set_isolated_buses_oos(bus_in_service, ppc_bus_isolated, net["_pd2ppc_lookups"]["bus"])
    #    mode = net["_options"]["mode"]
    elements = ["load", "motor", "sgen", "asymmetric_load", "asymmetric_sgen", "gen",
                "ward", "xward", "shunt", "ext_grid", "storage", "svc", "ssc"]  # ,"impedance_load"
    is_elements = dict()
    for element in elements:
        len_ = len(net[element].index)
        element_in_service = np.zeros(len_, dtype=bool)
        if len_ > 0:
            element_df = net[element]
            set_elements_oos(element_df["bus"].values, element_df["in_service"].values,
                             bus_in_service, element_in_service)
        if net["_options"]["mode"] == "opf" and element in ["load", "sgen", "storage"]:
            if "controllable" in net[element]:
                controllable = net[element].controllable.fillna(False).values.astype(bool)
                controllable_is = controllable & element_in_service
                if controllable_is.any():
                    is_elements["%s_controllable" % element] = controllable_is
                    element_in_service = element_in_service & ~controllable_is
        is_elements[element] = element_in_service

    is_elements["bus_is_idx"] = net["bus"].index.values[bus_in_service[net["bus"].index.values]]
    is_elements["line_is_idx"] = net["line"].index[net["line"].in_service.values]
    return is_elements
