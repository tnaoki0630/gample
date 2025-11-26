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
    # processes
    proc = []
    proc.append("artificialIonization")
    proc.append("hollowCathode")
    proc.append("integCDens_electron")
    proc.append("integCDens_ion_Xe1")
    proc.append("reduce_electron")
    proc.append("reduce_ion_Xe1")
    proc.append("solvePoisson")
    proc.append("update_electron")
    proc.append("update_ion_Xe1")
    # get array
    # logfiles = ['fine_20250801025821_log.xml']
    # logfiles = ['fine_3rd_20250902130112_log.xml']
    logfiles = ['org_20251004154620_log.xml']
    cycle_arr = []
    vals = []
    vals_arr = [[] for _ in range(9)]
    for file in logfiles:
        for i in range(9):
            cycles, val = parseXML(file, 'elapsedTime', proc[i])
            vals.append(val)
        cycle_arr.extend(cycles)
        for i, val in enumerate(vals):
            vals_arr[i].extend(val)
    # plot array
    values = vals_arr
    labels = proc
    weight = cycle_arr[1]-cycle_arr[0]
    scales = [1 for _ in range(9)]
    plot_values(cycle_arr, values, labels, scales, "elapsed time [s]", dt=5e-6)
    total = 0
    for group in zip(*vals_arr):
        total += sum(group)
    print("dt = ", weight*5e-3, " [ns], cycles = ", len(cycle_arr)*weight, " [cycles], elapsed time = ", total*weight*1e-6, " [sec], speed = ", total*1e-3/len(cycle_arr), " [msec/cyc]")
