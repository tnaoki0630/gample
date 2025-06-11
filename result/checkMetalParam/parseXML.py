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
    cycles, tgs1_ics128 = parseXML('tgs1ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs2_ics128 = parseXML('tgs2ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs4_ics128 = parseXML('tgs4ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs8_ics128 = parseXML('tgs8ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs16_ics128 = parseXML('tgs16ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs32_ics128 = parseXML('tgs32ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs64_ics128 = parseXML('tgs64ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs128_ics128 = parseXML('tgs128ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    timeICDE = [ tgs1_ics128 ,  tgs2_ics128 ,  tgs4_ics128 ,  tgs8_ics128 ,  tgs16_ics128 ,  tgs32_ics128 ,  tgs32_ics128 ,  tgs64_ics128 ,  tgs128_ics128 ]
    labels =   ["tgs1_ics128", "tgs2_ics128", "tgs4_ics128", "tgs8_ics128", "tgs16_ics128", "tgs32_ics128", "tgs32_ics128", "tgs64_ics128", "tgs128_ics128"]
    plot_values(cycles, timeICDE, labels, "elapsed time [msec]", 1e-3)
    
    cycles, tgs1_ics1 = parseXML('tgs1ics1_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs1_ics2 = parseXML('tgs1ics2_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs1_ics4 = parseXML('tgs1ics4_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs1_ics8 = parseXML('tgs1ics8_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs1_ics16 = parseXML('tgs1ics16_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs1_ics32 = parseXML('tgs1ics32_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs1_ics64 = parseXML('tgs1ics64_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs1_ics128 = parseXML('tgs1ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs1_ics256 = parseXML('tgs1ics256_log.xml', 'elapsedTime', 'integCDens_electron')
    timeICDE = [ tgs1_ics4 , tgs1_ics8 , tgs1_ics16 , tgs1_ics32 ,  tgs1_ics64 ,  tgs1_ics128 ,  tgs1_ics256 ]
    labels =   ["tgs1_ics4","tgs1_ics8","tgs1_ics16","tgs1_ics32", "tgs1_ics64", "tgs1_ics128", "tgs1_ics256"]
    plot_values(cycles, timeICDE, labels, "elapsed time [msec]", 1e-3)
    
    cycles, tgs2_ics2 = parseXML('tgs2ics2_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs2_ics4 = parseXML('tgs2ics4_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs2_ics8 = parseXML('tgs2ics8_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs2_ics16 = parseXML('tgs2ics16_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs2_ics32 = parseXML('tgs2ics32_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs2_ics64 = parseXML('tgs2ics64_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs2_ics128 = parseXML('tgs2ics128_log.xml', 'elapsedTime', 'integCDens_electron')
    cycles, tgs2_ics256 = parseXML('tgs2ics256_log.xml', 'elapsedTime', 'integCDens_electron')
    timeICDE = [ tgs2_ics4 , tgs2_ics8 , tgs2_ics16 , tgs2_ics32 ,  tgs2_ics64 ,  tgs2_ics128 ,  tgs2_ics256 ]
    labels =   ["tgs2_ics4","tgs2_ics8","tgs2_ics16","tgs2_ics32", "tgs2_ics64", "tgs2_ics128", "tgs2_ics256"]
    plot_values(cycles, timeICDE, labels, "elapsed time [msec]", 1e-3)
    
    cycles, tgs1_ics128 = parseXML('tgs1ics128_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs2_ics128 = parseXML('tgs2ics128_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs4_ics128 = parseXML('tgs4ics128_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs8_ics128 = parseXML('tgs8ics128_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs16_ics128 = parseXML('tgs16ics128_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs32_ics128 = parseXML('tgs32ics128_log.xml', 'memoryUsage', 'physicalFootprint')
    timeICDE = [ tgs1_ics128 ,  tgs2_ics128 ,  tgs4_ics128 ,  tgs8_ics128 ,  tgs16_ics128 ,  tgs32_ics128 ]
    labels =   ["tgs1_ics128", "tgs2_ics128", "tgs4_ics128", "tgs8_ics128", "tgs16_ics128", "tgs32_ics128"]
    plot_values(cycles, timeICDE, labels, "memory usage [MB]", 1e-3)

    cycles, tgs1_ics32 = parseXML('tgs1ics32_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs1_ics64 = parseXML('tgs1ics64_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs1_ics128 = parseXML('tgs1ics128_log.xml', 'memoryUsage', 'physicalFootprint')
    cycles, tgs1_ics256 = parseXML('tgs1ics256_log.xml', 'memoryUsage', 'physicalFootprint')
    timeICDE = [ tgs1_ics32 ,  tgs1_ics64 ,  tgs1_ics128 ,  tgs1_ics256 ]
    labels =   ["tgs1_ics32", "tgs1_ics64", "tgs1_ics128", "tgs1_ics256"]
    plot_values(cycles, timeICDE, labels, "memory usage [MB]", 1e-3)

