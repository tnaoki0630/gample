import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import sys

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
            # else:
            #     print(section_name,sec.get('Name'))
                

        if section is None:
            # No section found for this cycle
            values.append(None)
            continue

        # Single-tag case
        elem = section.find(tag_name)
        val = float(elem.text) if elem is not None and elem.text else None
        values.append(val)

        # if(tag_name == "totalEEincrement"):wEEincrement"),elem.text)

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
    # # processes
    # proc = []
    # proc.append("artificialIonization")
    # proc.append("hollowCathode")
    # proc.append("integCDens_electron")
    # proc.append("integCDens_ion_Xe1")
    # proc.append("reduce_electron")
    # proc.append("reduce_ion_Xe1")
    # proc.append("solvePoisson")
    # proc.append("update_electron")
    # proc.append("update_ion_Xe1")

    args = sys.argv
    pname = args[1]
    
    # get array
    # logfiles = ['fine_20250801025821_log.xml']
    # logfiles = ['fine_3rd_20250902130112_log.xml']
    logfiles = [pname]
    cycle_arr = []
    totEE_arr = []
    totEEinc_arr = []
    totJH_arr = []
    totKE_ele_arr = []
    totKE_ion_arr = []
    totKEinc_ele_arr = []
    totKEinc_ion_arr = []
    for file in logfiles:
        cycles, totEE = parseXML(file, 'solvePoisson', "totalElectricEnergy")
        cycles, totEEinc = parseXML(file, 'solvePoisson', "totalEEincrement")
        cycles, totJH = parseXML(file, 'solvePoisson', "jouleHeating")
        cycles, totKE_ele = parseXML(file, 'flowout_electron', "totalKineticEnergy")
        cycles, totKEinc_ele = parseXML(file, 'flowout_electron', "totalKEincrement")
        cycles, totKE_ion = parseXML(file, 'flowout_ion_Xe1', "totalKineticEnergy")
        cycles, totKEinc_ion = parseXML(file, 'flowout_ion_Xe1', "totalKEincrement")
        # cycle_arr.extend(cycles)
        # totEE_arr.extend(totEE)
        # totEEinc_arr.extend(totEEinc)
        # totJH_arr.extend(totJH)
        # totKE_ele_arr.extend(totKE_ele)
        # totKEinc_ele_arr.extend(totKEinc_ele)
        # totKE_ion_arr.extend(totKE_ion)
        # totKEinc_ion_arr.extend(totKEinc_ion)
    # plot array
    values = [totEE, totKE_ele, totKE_ion]
    labels = ["totEE", "totKE_ele", "totKE_ion"]
    scales = [1e-9 for _ in range(len(values))]
    # weight = 1.250000e+03
    # scales = [1,1,weight,weight]
    plot_values(cycles, values, labels, scales, "energy [GeV]", dt=5e-6)
    
    # increment
    values = [totEEinc, totKEinc_ele, totKEinc_ion, totJH]
    labels = ["totEEinc", "totKEinc_ele", "totKEinc_ion", "totJH"]
    scales = [1e-6 for _ in range(len(values))]
    # weight = 1.250000e+03
    # scales = [1,1,weight,weight]
    plot_values(cycles, values, labels, scales, "energy increment [MeV]", dt=5e-6)
    
    ## check
    # KE_old = 0
    # for KE, KEinc in zip(totKE_ele,totKEinc_ele):
    #     print(KE,float(KE)-KE_old,KEinc)
    #     KE_old = float(KE)
