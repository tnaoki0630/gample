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
    plot_variable_from_xml('log_cp.xml', 'elapsedTime', 'all')
    plot_variable_from_xml('log_cp.xml', 'flowout_ion_Xe1', 'all')
    plot_variable_from_xml('log_cp.xml', 'flowout_electron', 'all')
    plot_variable_from_xml('log_cp.xml', 'solvePoisson', 'meanCathode')
    plot_variable_from_xml('log_cp.xml', 'solvePoisson', 'iteration')
    plot_variable_from_xml('log_cp.xml', 'memoryUsage', 'all')
    # plot_variable_from_xml('log_cp.xml', 'outputField', 'all')

    cycles, Ne = parseXML('log_cp.xml', 'flowout_electron', 'particleNumber')
    cycles, Ni = parseXML('log_cp.xml', 'flowout_ion_Xe1', 'particleNumber')
    values = [Ne,Ni]
    labels = ["Ne","Ni"]
    weight = 6.25e4 # 1/cm
    domainAR = (2.5*1.25)/(2.5*1.0) # mthesis vs benchmark
    scales = [weight*domainAR,weight*domainAR]
    plot_values(cycles, values, labels, scales, "particle Number [1/cm]", dt=5e-6, xmin=0, xmax=20, ymin=0, ymax=8e11)

    cycles, Gamma_ic = parseXML('log_cp.xml', 'flowout_ion_Xe1', 'Xmin')
    cycles, Gamma_ia = parseXML('log_cp.xml', 'flowout_ion_Xe1', 'Xmax')
    cycles, Gamma_ec = parseXML('log_cp.xml', 'flowout_electron', 'Xmin')
    cycles, Gamma_ea = parseXML('log_cp.xml', 'flowout_electron', 'Xmax')
    Ly = 0.005*200 # cm
    dt = 5e-12  # s
    Gamma_i = [v1+v2 for v1,v2 in zip(Gamma_ic, Gamma_ia)]
    Gamma_e = [v1+v2 for v1,v2 in zip(Gamma_ec, Gamma_ea)]
    values = [Gamma_e,Gamma_i]
    labels = ["Gamma_e","Gamma_i"]
    scales = [weight/Ly/dt,weight/Ly/dt]
    plot_values(cycles, values, labels, scales, "particle Flux [1/cm2s]", dt=5e-6, xmin=0, xmax=20, ymin=0, ymax=4e17)
    
    cycles, addn_HC = parseXML('log_cp.xml', 'hollowCathode', 'addn_injected')
    cycles, addn_HC_neg = parseXML('log_cp.xml', 'hollowCathode', 'addn_neglected')
    cycles, addn_AI = parseXML('log_cp.xml', 'artificialIonization', 'addn_injected')
    values = [addn_HC,addn_HC_neg,addn_AI]
    labels = ["addn_HC","addn_HC_neg","addn_AI"]
    scales = [1,1,1]
    plot_values(cycles, values, labels, scales, "injected particles [-]", dt=5e-6, xmin=0)