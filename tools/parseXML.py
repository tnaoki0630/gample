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
    plt.xlim(left=0)
    # plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    cycle = 200
    cycles, arr1 = parseXML('large4_light_woFC_log.xml', 'elapsedTime', 'update_electron')
    cycles, arr2 = parseXML('large4_light_woFC_log.xml', 'elapsedTime', 'update_ion_Xe1')
    cycles, arr3 = parseXML('large4_light_woFC_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, arr4 = parseXML('large4_light_woFC_log.xml', 'elapsedTime', 'integCDens_ion_Xe1')
    update_woFC = [v1+v2+v3+v4 for v1, v2, v3, v4 in zip(arr1,arr2,arr3,arr4)]
    cycles, arr1 = parseXML('large4_light_wFC_log.xml', 'elapsedTime', 'update_electron')
    cycles, arr2 = parseXML('large4_light_wFC_log.xml', 'elapsedTime', 'update_ion_Xe1')
    cycles, arr3 = parseXML('large4_light_wFC_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, arr4 = parseXML('large4_light_wFC_log.xml', 'elapsedTime', 'integCDens_ion_Xe1')
    update_wFC = [v1+v2+v3+v4 for v1, v2, v3, v4 in zip(arr1,arr2,arr3,arr4)]
    cycles, arr1 = parseXML('large4_light_margedKarnel_log.xml', 'elapsedTime', 'update_electron')
    cycles, arr2 = parseXML('large4_light_margedKarnel_log.xml', 'elapsedTime', 'update_ion_Xe1')
    update_margedKarnel = [v1+v2 for v1, v2 in zip(arr1,arr2)]
    values = [update_woFC, update_wFC, update_margedKarnel]
    labels = ["updateAndDeposition_woFC", "updateAndDeposition_wFC", "update_margedKarnel"]
    scales = [1e-3,1e-3,1e-3]
    plot_values(cycles[0:cycle], values, labels, scales, "elapsed time [ms]", dt=5e-3)