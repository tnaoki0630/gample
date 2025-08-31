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

def plot_values(cycles, arr_values, labels, scales, title="", dt=1.0, xmin=None, xmax=None, ymin=None, ymax=None):
    plt.figure()
    for values, label_values, scale in zip(arr_values, labels, scales):
        plt.plot([cyc*dt for cyc in cycles], [val*scale for val in values], label=label_values)
    plt.legend(loc='lower right')
    plt.margins(x=0.0, y=0.30)
    plt.ylabel(title)
    plt.xlabel('time [us]')
    if xmin is not None:
        plt.xlim(left=xmin)
    if xmax is not None:
        plt.xlim(right=xmax)
    if ymin is not None:
        plt.ylim(bottom=ymin)
    if ymax is not None:
        plt.ylim(top=ymax)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # get array
    cycle_arr = []
    Ne_arr = []
    Ni_arr = []
    Gammae_arr = []
    Gammai_arr = []
    logfiles = ['fine_20250730223753_log.xml']
    logfiles.append('fine_20250801025821_log.xml')
    logfiles.append('fine_20250802040148_log.xml')
    logfiles.append('log_cp.xml')
    for file in logfiles:
        cycles, Ne = parseXML(file, 'flowout_electron', 'particleNumber')
        cycles, Ni = parseXML(file, 'flowout_ion_Xe1', 'particleNumber')
        cycles, Gammae = parseXML(file, 'flowout_electron', 'Xmax')
        cycles, Gammai = parseXML(file, 'flowout_ion_Xe1', 'Xmax')
        cycle_arr.extend(cycles)
        Ne_arr.extend(Ne)
        Ni_arr.extend(Ni)
        Gammae_arr.extend(Gammae)
        Gammai_arr.extend(Gammai)
    # plot array
    values = [Ne_arr, Ni_arr]
    labels = ["Ne", "Ni"]
    weight = 7.812500e3
    scales = [weight*1.0/1.25/4, weight*1.0/1.25/4]
    plot_values(cycle_arr, values, labels, scales, "particle Number [1/cm]", dt=5e-6, xmin=0, xmax=20, ymin=0)
    Ly = 0.005*200 # cm
    dt = 5e-12  # s
    values = [[-v for v in Gammae_arr],Gammai_arr]
    labels = ["Gamma_e","Gamma_i"]
    scales = [weight/Ly/dt,weight/Ly/dt]
    plot_values(cycle_arr, values, labels, scales, "particle Flux [1/cm2s]", dt=5e-6, xmin=0, xmax=20, ymin=0, ymax=4e17)