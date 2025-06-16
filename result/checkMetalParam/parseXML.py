import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

def parseXML(xml_file, section_name, tag_name):
    """
    Parses the given XML file.
    usage:: obtain time-series of elapsed time in integrate charge-density for electron.
    cycles, values = parseXML('log_log.xml', 'elapsedTime', 'integCDens_electron')
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
    plt.ylabel(title)
    plt.xlabel('Cycle')
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()    

if __name__ == '__main__':
    cycles, arr1 = parseXML('tgs1_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, arr2 = parseXML('tgs1_log.xml', 'elapsedTime', 'integCDens_ion_Xe1')
    cycles, arr3 = parseXML('tgs1_log.xml', 'elapsedTime', 'update_electron')
    cycles, arr4 = parseXML('tgs1_log.xml', 'elapsedTime', 'update_ion_Xe1')
    tgs1 = [v1+v2+v3+v4 for v1, v2, v3, v4 in zip(arr1,arr2,arr3,arr4)]
    cycles, arr1 = parseXML('tgs4_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, arr2 = parseXML('tgs4_log.xml', 'elapsedTime', 'integCDens_ion_Xe1')
    cycles, arr3 = parseXML('tgs4_log.xml', 'elapsedTime', 'update_electron')
    cycles, arr4 = parseXML('tgs4_log.xml', 'elapsedTime', 'update_ion_Xe1')
    tgs4 = [v1+v2+v3+v4 for v1, v2, v3, v4 in zip(arr1,arr2,arr3,arr4)]
    cycles, arr1 = parseXML('tgs16_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, arr2 = parseXML('tgs16_log.xml', 'elapsedTime', 'integCDens_ion_Xe1')
    cycles, arr3 = parseXML('tgs16_log.xml', 'elapsedTime', 'update_electron')
    cycles, arr4 = parseXML('tgs16_log.xml', 'elapsedTime', 'update_ion_Xe1')
    tgs16 = [v1+v2+v3+v4 for v1, v2, v3, v4 in zip(arr1,arr2,arr3,arr4)]
    cycles, arr1 = parseXML('tgs64_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, arr2 = parseXML('tgs64_log.xml', 'elapsedTime', 'integCDens_ion_Xe1')
    cycles, arr3 = parseXML('tgs64_log.xml', 'elapsedTime', 'update_electron')
    cycles, arr4 = parseXML('tgs64_log.xml', 'elapsedTime', 'update_ion_Xe1')
    tgs64 = [v1+v2+v3+v4 for v1, v2, v3, v4 in zip(arr1,arr2,arr3,arr4)]
    cycles, arr1 = parseXML('tgs128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, arr2 = parseXML('tgs128_log.xml', 'elapsedTime', 'integCDens_ion_Xe1')
    cycles, arr3 = parseXML('tgs128_log.xml', 'elapsedTime', 'update_electron')
    cycles, arr4 = parseXML('tgs128_log.xml', 'elapsedTime', 'update_ion_Xe1')
    tgs128 = [v1+v2+v3+v4 for v1, v2, v3, v4 in zip(arr1,arr2,arr3,arr4)]
    cycles, arr1 = parseXML('tgs256_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, arr2 = parseXML('tgs256_log.xml', 'elapsedTime', 'integCDens_ion_Xe1')
    cycles, arr3 = parseXML('tgs256_log.xml', 'elapsedTime', 'update_electron')
    cycles, arr4 = parseXML('tgs256_log.xml', 'elapsedTime', 'update_ion_Xe1')
    tgs256 = [v1+v2+v3+v4 for v1, v2, v3, v4 in zip(arr1,arr2,arr3,arr4)]
    timeICDE = [ tgs1 ,  tgs4 ,  tgs16 ,  tgs64 ,  tgs128 ,  tgs256 ]
    labels =   ["tgs1", "tgs4", "tgs16", "tgs64", "tgs128", "tgs256"]
    plot_values(cycles, timeICDE, labels, "elapsed time [msec]", 1e-3)
    
    cycles, tgs1 = parseXML('tgs1_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs2 = parseXML('tgs2_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs4 = parseXML('tgs4_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs8 = parseXML('tgs8_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs16 = parseXML('tgs16_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs32 = parseXML('tgs32_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs64 = parseXML('tgs64_log.xml', 'memoryUsage', 'physicalFootprint')
    timeICDE = [ tgs1 ,  tgs2 ,  tgs4 ,  tgs8 ,  tgs16 ,  tgs32 ,  tgs64 ]
    labels =   ["tgs1", "tgs2", "tgs4", "tgs8", "tgs16", "tgs32", "tgs64"]
    plot_values(cycles, timeICDE, labels, "memory usage [MB]", 1e-3)


