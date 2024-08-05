import pandapower as pp
import pandapower.shortcircuit as sc

def three_bus_example():
    net = pp.create_empty_network(sn_mva=100)
    b1 = pp.create_bus(net, 110)
    b2 = pp.create_bus(net, 110)
    b3 = pp.create_bus(net, 110)
    #b4 = pp.create_bus(net, 110)

    pp.create_ext_grid(net, b1, vm_pu=1.0,  s_sc_max_mva=100., s_sc_min_mva=80., rx_min=0.4, rx_max=0.4)
    pp.create_line(net, b1, b2, std_type="305-AL1/39-ST1A 110.0" , length_km=20.)
    pp.create_line(net, b2, b3, std_type="N2XS(FL)2Y 1x185 RM/35 64/110 kV" , length_km=15.)
    #pp.create_line(net, b3, b4, std_type="N2XS(FL)2Y 1x185 RM/35 64/110 kV", length_km=15.)
    #net.line["endtemp_degree"] = 80

    #pp.create_load(net, bus=b3, p_mw=5, q_mvar=1)

    #pp.create_sgen(net, b2, sn_mva=2, p_mw=0, k=1.2)

    net.ext_grid['x0x_min'] = 0.1
    net.ext_grid['r0x0_min'] = 0.1
    net.ext_grid['x0x_max'] = 0.1
    net.ext_grid['r0x0_max'] = 0.1

    net.line['r0_ohm_per_km'] = 0.1
    net.line['x0_ohm_per_km'] = 0.1
    net.line['c0_nf_per_km'] = 0.1
    net.line["endtemp_degree"] = 80
    return net

net = three_bus_example()

sc.calc_sc(net, case="max", bus=1, fault='1ph', branch_results=True, return_all_currents=True)