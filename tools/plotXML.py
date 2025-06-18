import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

def plot_variable_from_xml(xml_file, section_name, tag_name, ):
    """
    Parses the given XML file.
    """
    # Parse XML
    tree = ET.parse(xml_file)
    root = tree.getroot()

    cycle_ids = []
    # Containers for data
    if tag_name.lower() == 'all':
        tag_names = None
        tag_values = {}
    else:
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
            if tag_name.lower() == 'all' and tag_names:
                for tn in tag_names:
                    tag_values[tn].append(None)
            elif tag_name.lower() != 'all':
                values.append(None)
            continue

        if tag_name.lower() == 'all':
            # Initialize tag_names and lists on first encounter
            if tag_names is None:
                tag_names = [elem.tag for elem in section]
                for tn in tag_names:
                    tag_values[tn] = []
            # Append values for each tag
            for tn in tag_names:
                elem = section.find(tn)
                val = float(elem.text) if elem is not None and elem.text else None
                tag_values[tn].append(val)
        else:
            # Single-tag case
            elem = section.find(tag_name)
            val = float(elem.text) if elem is not None and elem.text else None
            values.append(val)

    # Plotting
    plt.figure()
    if tag_name.lower() == 'all':
        for tn in tag_names:
            # 粒子数は桁が違いすぎるので除外
            if tn != "particleNumber": plt.plot(cycle_ids, tag_values[tn], label=tn)
        plt.ylabel('Value')
    else:
        plt.plot(cycle_ids, values, label=tag_name)
        plt.ylabel(tag_name)
    plt.legend(loc='lower right')
    plt.xlabel('Cycle ID')
    plt.title(f'{section_name} : {tag_name} vs Cycle ID')
    plt.tight_layout()
    plt.show()

# Example usage:
# plot_variable_from_xml('log.xml', 'solvePoisson', 'iteration')
# plot_variable_from_xml('log.xml', 'elapsedTime', 'all')

if __name__ == '__main__':
    plot_variable_from_xml('log_cp.xml', 'elapsedTime', 'all')
    plot_variable_from_xml('log_cp.xml', 'injection_electron', 'all')
    plot_variable_from_xml('log_cp.xml', 'flowout_ion_Xe1', 'all')
    plot_variable_from_xml('log_cp.xml', 'flowout_electron', 'all')
    plot_variable_from_xml('log_cp.xml', 'flowout_ion_Xe1', 'particleNumber')
    plot_variable_from_xml('log_cp.xml', 'flowout_electron', 'particleNumber')
    plot_variable_from_xml('log_cp.xml', 'solvePoisson', 'meanCathode')
    plot_variable_from_xml('log_cp.xml', 'solvePoisson', 'iteration')
    plot_variable_from_xml('log_cp.xml', 'memoryUsage', 'all')
    plot_variable_from_xml('log_cp.xml', 'outputField', 'all')
