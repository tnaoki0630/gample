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
    fig, ax = plt.subplots()
    for values, label_values in zip(arr_values, labels):
        ax.plot(cycles, [val*scale for val in values], label=label_values)
    ax.margins(x=0.0, y=0.30)
    ax.legend(loc='lower right')
    plt.ylabel(title)
    plt.xlabel('Cycle')
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    cycles, timeICDE = parseXML('log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, timeICDE_before = parseXML('log_before.xml', 'elapsedTime', 'integCDens_electron')
    timeICDE = [ timeICDE , timeICDE_before ]
    labels =   ["timeICDE","timeICDE_before"]
    plot_values(cycles, timeICDE, labels, "elapsed time [msec]", 1e-3)
    cycles, physUsage = parseXML('log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, physUsage_before = parseXML('log_before.xml', 'memoryUsage', 'physicalFootprint')
    timeICDE = [ physUsage , physUsage_before ]
    labels =   ["physUsage","physUsage_before"]
    plot_values(cycles, timeICDE, labels, "physical memory usage [MB]", 1e-3)