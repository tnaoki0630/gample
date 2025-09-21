from pathlib import Path
import xml.etree.ElementTree as ET

xml_dir = Path("docs/doxygen/xml")
for f in xml_dir.glob("*.xml"):
    tree = ET.parse(f)
    root = tree.getroot()

    def prune(e):
        for child in list(e):
            if child.tag == "memberdef" and child.attrib.get("kind") == "property":
                e.remove(child)
            else:
                prune(child)
    prune(root)
    tree.write(f, encoding="utf-8", xml_declaration=True)
