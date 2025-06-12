import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

def parseXML(xml_file, section_name, tag_name):
    """
    Parses the given XML file.
    usage:: obtain time-series of elapsed time in integrate charge-density for electron.
    cycles, values = parseXML('log.xml', 'elapsedTime', 'integCDens_electron')
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

def plot_values(cycles, arr_values, labels, title="", scale=1.0):
    plt.figure()
    for values, label_values in zip(arr_values, labels):
        plt.plot(cycles, [val*scale for val in values], label=label_values)
    plt.legend(loc='lower right')
    plt.margins(x=0.0, y=0.30)
    plt.ylabel(title)
    plt.xlabel('Cycle')
    # plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    cycles, OMP1 = parseXML('CG_OMP1_log.xml', 'elapsedTime', 'solvePoisson')
    cycles, OMP2 = parseXML('CG_OMP2_log.xml', 'elapsedTime', 'solvePoisson')
    cycles, OMP3 = parseXML('CG_OMP3_log.xml', 'elapsedTime', 'solvePoisson')
    cycles, OMP4 = parseXML('CG_OMP4_log.xml', 'elapsedTime', 'solvePoisson')
    cycles, OMP5 = parseXML('CG_OMP5_log.xml', 'elapsedTime', 'solvePoisson')
    cycles, OMP6 = parseXML('CG_OMP6_log.xml', 'elapsedTime', 'solvePoisson')
    values = [OMP1, OMP2, OMP3, OMP4, OMP5, OMP6]
    labels = ["CG_OMP1", "CG_OMP2", "CG_OMP3", "CG_OMP4", "CG_OMP5", "CG_OMP6"]
    plot_values(cycles, values, labels, "elapsed time [ms]", 1e-3)
    
    cycles, OMP1 = parseXML('CG_OMP1_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, OMP2 = parseXML('CG_OMP2_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, OMP3 = parseXML('CG_OMP3_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, OMP4 = parseXML('CG_OMP4_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, OMP5 = parseXML('CG_OMP5_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, OMP6 = parseXML('CG_OMP6_log.xml', 'elapsedTime', 'integCDens_electron')
    values = [OMP1, OMP2, OMP3, OMP4, OMP5, OMP6]
    labels = ["CG_OMP1", "CG_OMP2", "CG_OMP3", "CG_OMP4", "CG_OMP5", "CG_OMP6"]
    plot_values(cycles, values, labels, "elapsed time [ms]", 1e-3)
    
    # cycles, CG = parseXML('CG_log.xml', 'solvePoisson', 'iteration')
    # cycles, BiCGStab = parseXML('BiCGStab_log.xml', 'solvePoisson', 'iteration')
    # values = [CG, BiCGStab]
    # labels = ["CG", "BiCGStab"]
    # plot_values(cycles, values, labels, "iterations", 1.0)