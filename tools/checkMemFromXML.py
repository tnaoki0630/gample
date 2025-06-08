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

def plot_2values(
    cycles,
    arr_values=None,
    labels=None,
    title="",
    scale=1.0,
    arr_values2=None,
    labels2=None,
    title2="",
    scale2=1.0
):
    """
    cycles:    x 軸に対応するリスト（一例: [0,1,2,...]）
    arr_values:   左軸にプロットしたい数値系列のリスト（例: [[v11, v12,...], [v21, v22,...], ...]）
    labels:    arr_values に対応する凡例ラベルのリスト
    title:     左軸のラベル（Y 軸タイトル）
    scale:     arr_values に乗ずるスケーリング係数
    arr_values2:  右軸にプロットしたい数値系列のリスト
    labels2:   arr_values2 に対応する凡例ラベルのリスト
    title2:    右軸のラベル
    scale2:    arr_values2 に乗ずるスケーリング係数
    """
    if arr_values is None:
        arr_values = []
    if labels is None:
        labels = []
    if arr_values2 is None:
        arr_values2 = []
    if labels2 is None:
        labels2 = []

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    cmap = plt.get_cmap("tab10")

    # 左軸のプロット
    for i, vals in enumerate(arr_values):
        # vals: 1系列分の数値リスト
        scaled = [v * scale for v in vals]
        ax1.plot(cycles, scaled, label=labels[i] if i < len(labels) else None, color=cmap(i))

    # 右軸のプロット
    offset = len(arr_values)
    for j, vals2 in enumerate(arr_values2):
        scaled2 = [v * scale2 for v in vals2]
        ax2.plot(cycles, scaled2, label=labels2[j] if j < len(labels2) else None, color=cmap(j + offset), linestyle='--')

    # 凡例を両方の軸からまとめて取得
    lines1, lab1 = ax1.get_legend_handles_labels()
    lines2, lab2 = ax2.get_legend_handles_labels()
    all_lines = lines1 + lines2
    all_labels = lab1 + lab2
    if all_lines:
        ax1.legend(all_lines, all_labels, loc='lower right')

    # 軸ラベル、タイトル、制限など
    ax1.set_ylabel(title)
    ax2.set_ylabel(title2)
    ax1.set_xlabel('Cycle')
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax2.set_ylim(bottom=0)
    plt.show()


if __name__ == '__main__':
    cycles, memUsage_physFP = parseXML('log_cp.xml', 'memoryUsage', 'physicalFootprint')
    cycles, memUsage_resiSS = parseXML('log_cp.xml', 'memoryUsage', 'residentSetSize')
    cycles, pNum_electron = parseXML('log_cp.xml', 'flowout_electron', 'particleNumber')
    cycles, pNum_ion_Xe1 = parseXML('log_cp.xml', 'flowout_ion_Xe1', 'particleNumber')
    values = [memUsage_physFP, memUsage_resiSS]
    values2 = [pNum_electron, pNum_ion_Xe1]
    labels =   ["memUsage_physFP", "memUsage_resiSS"]
    labels2 = ["pNum_electron", "pNum_ion_Xe1"]
    plot_2values(cycles, values, labels, "memory usage [MB]", 1e-3, values2, labels2, "particle number [1e6]", 1e-6)
