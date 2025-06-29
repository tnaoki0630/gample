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

def plot_values(cycles, arr_values, labels, scales, title="", xscale=1.0):
    plt.figure()
    for values, label_values, scale in zip(arr_values, labels, scales):
        plt.plot([cyc*xscale for cyc in cycles], [val*scale for val in values], label=label_values)
    plt.legend(loc='lower right')
    plt.margins(x=0.0, y=0.30)
    plt.ylabel(title)
    plt.xlabel('time [ns]')
    # plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    cycle = 3200
    cycles, Ne_long = parseXML('long_20250620033134_log.xml', 'flowout_electron', 'particleNumber')
    cycles, Ni_long = parseXML('long_20250620033134_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    cycles, Ne_large = parseXML('large_20250621024603_log.xml', 'flowout_electron', 'particleNumber')
    cycles, Ni_large = parseXML('large_20250621024603_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    cycles, Ne_large2 = parseXML('large2_20250623042317_log.xml', 'flowout_electron', 'particleNumber')
    cycles, Ni_large2 = parseXML('large2_20250623042317_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    cycles, Ne_large3 = parseXML('large3_20250624212336_log.xml', 'flowout_electron', 'particleNumber')
    cycles, Ni_large3 = parseXML('large3_20250624212336_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    cycles, Ne_large4 = parseXML('large4_20250626124908_log.xml', 'flowout_electron', 'particleNumber')
    cycles, Ni_large4 = parseXML('large4_20250626124908_log.xml', 'flowout_ion_Xe1', 'particleNumber')
    values = [Ne_long[0:cycle], Ni_long[0:cycle], Ne_large[0:cycle], Ni_large[0:cycle], Ne_large2[0:cycle], Ni_large2[0:cycle], Ne_large3[0:cycle], Ni_large3[0:cycle], Ne_large4[0:cycle], Ni_large4[0:cycle]]
    labels = ["Ne_ippc10", "Ni_ippc10", "Ne_ippc20", "Ni_ippc20", "Ne_ippc30", "Ni_ippc30", "Ne_ippc50", "Ni_ippc50", "Ne_ippc100(ngy=50)", "Ni_ippc100(ngy=50)"]
    scales = [12.5e4, 12.5e4, 6.25e4, 6.25e4, 4.166667e4, 4.166667e4, 2.5e4, 2.5e4, 1.25e4*4, 1.25e4*4]
    plot_values(cycles[0:cycle], values, labels, scales, "particle Number [1/cm]", 5e-3)