import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

def parseXML(xml_file, section_name, tag_name):
    """
    Parses the given XML file.
    usage:: obtain time-series of elapsed time in integrate charge-density for electron.
    cycles, values = parseXML('log.xml', 'flowout_electron', 'integCDens_electron')
    """
    # Parse XML
    tree = ET.parse(xml_file)
    root = tree.getroot()

    cycle_ids = []
    values = []

    # Iterate over each Cycle element
    for cycle in root.findall('Cycle'):
        cycle_id = int(cycle.get('ID'))
        cycle_ids.append(cycle_id)

        # Find the matching Section
        section = None
        for sec in cycle.findall('Section'):
            if sec.get('Name') == section_name:
                section = sec
                break

        if section is None:
            # No section found for this cycle
            values.append(None)
            continue

        # Single-tag case
        elem = section.find(tag_name)
        val = float(elem.text) if elem is not None and elem.text else None
        values.append(val)

    return cycle_ids, values

def plot_values(cycles, arr_values, labels, scales, title="", dt=1.0):
    plt.figure()
    for values, label_values, scale in zip(arr_values, labels, scales):
        plt.plot([cyc*dt for cyc in cycles], [val*scale for val in values], label=label_values)
    plt.legend(loc='lower right')
    plt.margins(x=0.0, y=0.30)
    plt.ylabel(title)
    plt.xlabel('time [ns]')
    # plt.xlim(left=0)
    # plt.ylim(bottom=0)q
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    cycle = 3200
    #cycles, Ne_long = parseXML('long_20250619025651_log.xml', 'flowout_electron', 'particleNumber')
    #cycles, Ni_long = parseXML('long_20250619025651_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    # cycles, Ne_large = parseXML('large_20250621024603_log.xml', 'flowout_electron', 'particleNumber')
    # cycles, Ni_large = parseXML('large_20250621024603_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    # cycles, Ne_large2 = parseXML('large2_20250623042317_log.xml', 'flowout_electron', 'particleNumber')
    # cycles, Ni_large2 = parseXML('large2_20250623042317_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    # cycles, Ne_large3 = parseXML('large3_20250624212336_log.xml', 'flowout_electron', 'particleNumber')
    # cycles, Ni_large3 = parseXML('large3_20250624212336_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    # cycles, Ne_large4 = parseXML('large4_20250626124908_log.xml', 'flowout_electron', 'particleNumber')
    # cycles, Ni_large4 = parseXML('large4_20250626124908_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    # cycles, Ne_large5 = parseXML('log_cp.xml', 'flowout_electron', 'particleNumber')
    # cycles, Ni_large5 = parseXML('log_cp.xml', 'flowout_ion_Xe1', 'particleNumber')
    # values = [Ne_large[0:cycle], Ni_large[0:cycle], Ne_large2[0:cycle], Ni_large2[0:cycle], Ne_large3[0:cycle], Ni_large3[0:cycle], Ne_large4[0:cycle], Ni_large4[0:cycle], Ne_large5[0:cycle], Ni_large5[0:cycle]]
    # labels = ["Ne_ippc20", "Ni_ippc20", "Ne_ippc30", "Ni_ippc30", "Ne_ippc50", "Ni_ippc50", "Ne_ippc100(ngy=50)", "Ni_ippc100(ngy=50)", "Ne_ippc250(ngy=100,1st)", "Ni_ippc250(ngy=100,1st)"]
    # scales = [6.25e4, 6.25e4, 4.166667e4, 4.166667e4, 2.5e4, 2.5e4, 1.25e4*4, 1.25e4*4, 0.5e4*2, 0.5e4*2]
    # plot_values(cycles[0:cycle], values, labels, scales, "particle Number [1/cm]", 5e-3)
    
    cycles, ele_Xmin = parseXML('log_cp.xml', 'flowout_electron', 'Xmin')
    cycles, ele_Xmax = parseXML('log_cp.xml', 'flowout_electron', 'Xmax')
    cycles, ele_pull = parseXML('log_cp.xml', 'flowout_electron', 'pulledPtclNum')
    cycles, ion_Xmin = parseXML('log_cp.xml', 'flowout_ion_Xe1', 'Xmin')
    cycles, ion_Xmax = parseXML('log_cp.xml', 'flowout_ion_Xe1', 'Xmax')
    cycles, ion_pull = parseXML('log_cp.xml', 'flowout_ion_Xe1', 'pulledPtclNum')
    cycles, injection = parseXML('log_cp.xml', 'injection_electron', 'hollow-cathode')
    injectionError = [max(0,-(v1+v2))-v3 for v1,v2,v3 in zip(ele_Xmin,ion_Xmin,injection)]
    pullError_ele = [v1+v2+v3 for v1,v2,v3 in zip(ele_Xmin,ele_Xmax,ele_pull)]
    pullError_ion = [v1+v2-v3 for v1,v2,v3 in zip(ion_Xmin,ion_Xmax,ion_pull)]
    values = [injectionError,pullError_ele,pullError_ion]
    labels = ["injectionError","pullError_ele","pullError_ion"]
    scales = [1,1,1]
    plot_values(cycles, values, labels, scales, "particle Number [-]", dt=5e-3)
    values = [ele_Xmin, ion_Xmin, injection]
    labels = ["ele_Xmin", "ion_Xmin", "injection"]
    scales = [1,1,1]
    plot_values(cycles, values, labels, scales, "particle Number [-]", dt=5e-3)
